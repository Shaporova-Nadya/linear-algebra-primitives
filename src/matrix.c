#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct CSR;

void csrVecMul(const CSR* matrix, const double* vector, double* result)
{
    for(int i = 0; i < matrix->nrows; i++){
        double sum = 0.0;
        int start = matrix->rowPtr[i];
        int end = matrix->rowPtr[i + 1];
        for(int j = start; j < end; j++){
            int col = matrix->colIndex[j];
            if(col < 0 || col >= matrix->ncols){
                printf("Ошибка: выход за диапазон\n");
                return -1;
            }
            sum += matrix->val[j] * vector[matrix->colIndex[j]];
        }
        result[i] = sum;
    }
}

CSR* csrMatMul(const CSR* A, const CSR* B)
{
    if(A->ncols != B->nrows){
        printf("Ошибка: матрицы нельзя перемножить");
        return NULL;
    }

    int m = A->nrows;
    int n = B->ncols;

    //Временные массивы для разреженного аккумулятора
    double* workspace = (double*)calloc(n, sizeof(double));//через calloc, т.к. workspace должен быть изначально заполнен нулями
    int* visited = (int*)malloc(n * sizeof(int));
    for(int i = 0; i < n; i++){
        visited[i] = -1;//заполняем -1, чтобы значение в visited отличалось от номера строки
    }

    //Динамические массивы для результата
    int cap = 1024;//В промышленных реализациях применяется редко из-за накладных расходов на realloc
    int* colInds = (int*)malloc(cap * sizeof(int));
    double* vals = (double*)malloc(cap * sizeof(double));
    int* rowPtr = (int*)malloc((m + 1) * sizeof(int));
    if(!workspace || !visited || !colInds || !vals || !rowPtr){
        //Обработка ошибок
        free(workspace);
        free(visited);
        free(colInds);
        free(vals);
        free(rowPtr);
        return NULL;
    }

    rowPtr[0] = 0;
    int nnz = 0;

    for(int i = 0; i < m; i++){
        int startA = A->rowPtr[i];
        int endA = A->rowPtr[i + 1];
        for(int ja = startA; ja < endA; ja++){
            double aVal = A->val[ja];
            int k = A->colIndex[ja]; //Номер строки в B

            int startB = B->rowPtr[k];
            int endB = B->rowPtr[k + 1];
            for(int jb = startB; jb < endB; jb++){
                double bVal = B->val[jb];
                int j = B->colIndex[jb];
                if(visited[j] != i){
                    visited[j] = i;
                    workspace[j] = aVal * bVal;
                } else {
                    workspace[j] += aVal * bVal;
                }
            }
        }

        for(int j = 0; j < n; j++){
            if(visited[j] == 1){
                if(nnz >= cap){
                    cap *= 2;
                    colInds = (int*)realloc(colInds, cap * sizeof(int));
                    vals = (double*)realloc(vals, cap * sizeof(double));
                    if(!colInds || !vals){
                        return NULL;
                    }
                    
                    colInds[nnz] = j;
                    vals[nnz] = workspace[j];
                    nnz++;
                }
            }
            rowPtr[i + 1] = nnz;
        }

        //Освобождение временных массивов
        free(workspace);
        free(visited);

        //Формирование структуры результата
        CSR* C = (CSR*)malloc(sizeof(CSR));
        if(!C){
            free(colInds);
            free(vals);
            free(rowPtr);
            return NULL;
        }
        C->nrows = m;
        C->ncols = n;
        C->nnz = nnz;
        C->rowPtr = rowPtr;
        C->val = vals;
        C->colIndex = colInds;

        return C;
    }
}

CSR* csrMap(const CSR* A, double (*f)(double))
{
    if(!A){
        return NULL;
    }

    CSR* B = (CSR*)malloc(sizeof(CSR));

    if(!B){
        return NULL;
    }

    B->nrows = A->nrows;
    B->ncols = A->ncols;

    int cap = 1024;
    int* colInds = (int*)malloc(cap * sizeof(int));
    double* val = (double*)malloc(cap * sizeof(double));
    int* rowPtr = (int*)malloc((A->nrows + 1) * sizeof(int));

    if(!colInds || !val || !rowPtr){
        free(colInds);
        free(val);
        free(rowPtr);
        free(B);
        return NULL;
    }

    int nnz = 0;//Счетчик не 0 элементов
    int rowPtr[0] = 0;

    for(int i = 0; i < A->nrows; i++){
        int start = A->rowPtr[i];
        int end = A->rowPtr[i + 1];

        for(int j = start; j < end; j++){
            double newVal = f(A->val[j]);
            if(newVal != 0.0){
                if(nnz > cap){
                    cap *= 2;
                    colInds = (int*)realloc(colInds, cap * sizeof(int));
                    val = (double*)realloc(val, cap * sizeof(double));
                    if(!colInds || !val){
                        free(colInds);
                        free(val);
                        free(rowPtr);
                        free(B);
                    }
                }
                colInds[nnz] = A->colIndex[j];//В массив номеров столбцов записывается номер столбца текущего элемента
                val[nnz] = newVal;//В массив значений на ту же позицию записывается значение
                nnz++;//увеличиваем счетчик
            }
        }
        rowPtr[i + 1] = nnz;
    }
    //Сжимаем массивы до точного размера
    colInds = (int*)realloc(colInds, nnz * sizeof(int));
    val = (double*)realloc(val, nnz * sizeof(double));

    B->nnz = nnz;
    B->val = val;
    B->colIndex = colInds;
    B->rowPtr = rowPtr;

    return B;
}
CSR* csrMap2(const CSR* A, const CSR* B, double(*g)(double, double))
{
    if(!B || !A || A->nrows != B->nrows || A->ncols != B->ncols){
        return NULL;
    }

    int rows = A->nrows;
    int cols = A->ncols;//Не используется напрямую, но для проверок

    //Динамические массивы для результа
    int cap = 1024;
    int* colInds = (int*)malloc(cap * sizeof(int));
    int* rowPtr = (int*)malloc((rows + 1) * sizeof(int));
    double* vals = (double*)malloc(cap * sizeof(double));

    if(!colInds || !vals || !rowPtr){
        free(colInds);
        free(rowPtr);
        free(vals);
        return NULL;
    }

    int nnz = 0;
    int rowPtr[0] = 0;

    for(int i = 0; i < rows; i++){
        //Получаем диапазон для строки i в A и B
        int aStart = A->rowPtr[i];
        int aEnd = A->rowPtr[i + 1];
        int bStart = B->rowPtr[i];
        int bEnd = B->rowPtr[i + 1];

        int aIdx = aStart;
        int bIdx = bStart;

        //Слияние двух отсортированных списков столбцов
        while(aIdx < aEnd && bIdx < bEnd){
            int colA = A->colIndex[aIdx];
            int colB = B->colIndex[bIdx];
            double val;//Для хранения значений вычесленного столбца
            int col;//ДЛя значений вычесленного столбца
            if(colA == colB){
                val = g(A->val[aIdx], B->val[bIdx]);
                col = colA;
                aIdx++;
                bIdx++;
            } else if(colA < colB){
                val = g(A->val[aIdx], 0.0);
                col = colA;
                aIdx++;
            } else {
                val = g(0.0, B->val[bIdx]);
                col = colB;
                bIdx++;
            }

            if(val != 0.0){
                if(nnz >= cap){
                    cap *= 2;
                    double* newVal = (double*)realloc(vals, cap * sizeof(double));
                    int* newCol = (int*)realloc(colInds, cap * sizeof(int));
                    if(!newCol || !newVal){
                        free(newCol ? newCol : colInds);
                        free(newVal ? newVal : vals);
                        free(rowPtr);

                        return NULL;
                    }
                    colInds = newCol;
                    vals = newVal;
                }
                colInds[nnz] = A->colIndex[aIdx];
                vals[nnz] = val;
                nnz++;
            }
            aIdx++;
        }

        while(bIdx < bEnd){
            double val = g(0.0, B->val[bIdx]);
            if(val != 0.0){
                if(nnz > cap){
                    cap *= 2;
                    double* newVal = (double*)realloc(vals, cap * sizeof(double));
                    int* newCol = (int*)realloc(colInds, cap * sizeof(int));
                    if(!newCol || !newVal){
                        free(newCol ? newCol : colInds);
                        free(newVal ? newVal : vals);
                        free(rowPtr);

                        return NULL;
                    }
                    colInds = newCol;
                    vals = newVal;
                }
                colInds[nnz] = A->colIndex[aIdx];
                vals[nnz] = val;
                nnz++;
            }
            bIdx++;
        }
        rowPtr[i + 1] = nnz;
    }

    int* finalCol = (int*)realloc(colInds, nnz * sizeof(int));
    double* finalVal = (double*)realloc(vals, nnz * sizeof(double));
    if(nnz > 0 && !finalCol || !finalVal){
        free(finalCol ? finalVal : colInds);
        free(finalVal ? finalVal : vals);
        free(rowPtr);
        return NULL;
    }

    if(finalCol){
        colInds = finalCol;
    }

    if(finalVal){
        vals = finalVal;
    }

    CSR* C = (CSR*)malloc(sizeof(CSR));
    if(!C){
        free(colInds);
        free(vals);
        free(rowPtr);

        return NULL;
    }
    C->colIndex = colInds;
    C->ncols = cols;
    C->nnz = nnz;
    C->nrows = cols;
    C->val = vals;
    C->rowPtr = rowPtr;

    return C;
}