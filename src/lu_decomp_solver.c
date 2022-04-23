#include "lu_decomp_solver.h"
#include <math.h>
#include <stdint.h>
#include <stddef.h>
#include <malloc.h>

#define MEPS_LU (1e-10)

int32_t LU_decomposition(int32_t dim, int32_t *pivot, double *matrix) {
    int32_t ret = (int32_t)LU_NG;
    int32_t row, col, k;
    int32_t ip, ip_tmp;
    double tmp_val, max_val;

    if ((NULL != pivot) && (NULL != matrix)) {
        for (k = 0; k < dim; k++) {
            pivot[k] = k;
        }

        for (k = 0; k < dim - 1; k++) {
            max_val = fabs(matrix[k * dim + k]);
            ip = k;

            for (row = k + 1; row < dim; row++) {
                tmp_val = fabs(matrix[row * dim + k]);

                if (max_val < tmp_val) {
                    max_val = tmp_val;
                    ip = row;
                }
            }

            // 対角成分の最大値が非常に小さい場合
            if (max_val < (double)MEPS_LU) {
                goto EXIT_LU_DECOMP;
            }
            else if (ip != k) {
                for (col = k; col < dim; col++) {
                    tmp_val = matrix[ip * dim + col];
                    matrix[ip * dim + col] = matrix[k * dim + col];
                    matrix[k * dim + col] = tmp_val;
                }
                ip_tmp = pivot[ip];
                pivot[ip] = pivot[k];
                pivot[k] = ip_tmp;

                for (col = 0; col < k; col++) {
                    tmp_val = matrix[k * dim + col];
                    matrix[k * dim + col] = matrix[ip * dim + col];
                    matrix[ip * dim + col] = tmp_val;
                }
            }

            for (row = k + 1; row < dim; row++) {
                tmp_val = matrix[row * dim + k] / matrix[k * dim + k];
                matrix[row * dim + k] = tmp_val;

                for (col = k + 1; col < dim; col++) {
                    matrix[row * dim + col] -= tmp_val * matrix[k * dim + col];
                }
            }
        }
        ret = (int32_t)LU_OK;
    }
EXIT_LU_DECOMP:

    return ret;
}

int32_t LU_solver(int32_t dim, int32_t *pivot, double *matrix, double *vec) {
    int32_t ret = (int32_t)LU_NG;
    int32_t row, col;
    double sum;
    double *y = NULL;

    if ((NULL != pivot) && (NULL != matrix) && (NULL != vec)) {
        y = (double *)malloc(sizeof(double) * dim);
        if (NULL == y) {
            goto EXIT_LU_SOLVER;
        }

        // 前進代入
        // Ly = vecからyを計算
        for (row = 0; row < dim; row++) {
            sum = vec[pivot[row]];

            for (col = 0; col < row; col++) {
                sum -= matrix[row * dim + col] * y[col];
            }
            y[row] = sum;
        }

        // 後退代入
        // Ux = yからxを計算
        for (row = dim - 1; row >= 0; row--) {
            sum = y[row];

            for (col = row + 1; col < dim; col++) {
                sum -= matrix[row * dim + col] * vec[col];
            }
            vec[row] = sum / matrix[row * dim + row];
        }
        ret = (int32_t)LU_OK;
    }
EXIT_LU_SOLVER:
    if (NULL != y) {
        free(y);
    }

    return ret;
}