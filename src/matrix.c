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