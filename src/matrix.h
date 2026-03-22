# pragma once

/**
 * @brief Разреженная матрица в формате CSR.
 * 
 * Хранит матрицу размером nrows x ncols с nnz ненулевыми элементами.
 * Массивы val и colIndex имеют длину nnz, rowPtr - длину nrows + 1.
 */
typedef struct CSR{
    int nrows;
    int ncols;
    int nnz;
    int* val;
    int* colIndex;
    int* rowPtr;
}CSR;

/**
 * @brief Умножение матрицы CSR на вектор.
 * 
 * @param[in]  matrix   Указатель на матрицу
 * @param[in]  vector   Входной вектор
 * @param[out] result   Выходной вектор
 */
void csrVecMul(const CSR* matrix, const double* vector, double* result);

/**
 * @brief Перемножение двух CSR-матриц 
 * 
 * @param[in] A  Первая матрица (размером m x p)
 * @param[in] B  Вторая матрица (размером p x n)
 * @result Указатель на новую матрицу C (размером m x n) или NULL при ошибке
 */
CSR* csrMatMul(const CSR* A, const CSR* B);