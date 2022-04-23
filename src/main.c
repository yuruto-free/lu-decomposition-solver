#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include "lu_decomp_solver.h"

#define RET_OK_MAIN (0)
#define RET_NG_MAIN (1)
#define MAX_TEST_DATA (2)
#define MAX_DIM (4)

/**
 * @brief 行列ベクトル積の計算
 * 
 * @param[in]  dim    次元数
 * @param[in]  matrix 行列
 * @param[in]  vec    ベクトル
 * @param[out] out    演算結果
 * @return RET_OK_MAIN 正常終了
 * @return RET_NG_MAIN 異常終了
 */
static int32_t mul_matvec(int32_t dim, double *matrix, double *vec, double *out);
/**
 * @brief 結果の出力
 * 
 * @param[in] dim 次元数
 * @param[in] idx パターンのインデックス
 * @param[in] out 出力データ
 * @return RET_OK_MAIN 正常終了
 * @return RET_NG_MAIN 異常終了
 */
static int32_t print_vec(int32_t dim, int32_t idx, double *out);

struct data_t {
    int32_t dim;
    int32_t pattern;
    double *matrix;
    double *exacts;
};

double matrix3x3[9] = {
    3.0, -2.0,  4.0,
    2.5,  1.5, -2.0,
    0.5, -5.0, -3.0,
};
double matrix4x4[16] = {
     0.5,  4.0, 6.0, -2.5,
     4.0, -1.0, 3.5, -2.0,
    -1.5,  0.5, 4.5, -0.75,
     1.5,  2.0, 5.0,  4.0,
};
double ans3[12] = { // dim: 3, pattern: 4
     3.0, 2.0, 1.0, // case 1
     1.0, 1.0, 1.5, // case 2
    -2.0, 1.0, 3.0, // case 3
     1.0, 2.0, 3.0, // case 4
};
double ans4[20] = { // dim: 4, pattern: 5
     1.0,  2.0,  3.0,  4.0, // case 1
    -1.0,  2.0, -3.0,  4.0, // case 2
     4.0, -2.0,  2.0, -4.0, // case 3
     3.0,  6.0,  1.0, -3.0, // case 4
     4.0,  3.0,  2.0,  1.0, // case 5
};

struct data_t TEST_DATA[MAX_TEST_DATA] = {
    {3, 4, matrix3x3, ans3},
    {4, 5, matrix4x4, ans4},
};

int main(int argc, char **argv) {
    int32_t ret;
    int32_t i, j;
    int32_t dim, pattern;
    struct data_t *target;
    int32_t pivot[MAX_DIM];
    double exact[MAX_DIM], vec[MAX_DIM], matrix[MAX_DIM*MAX_DIM];

    for (i = 0; i < (int32_t)MAX_TEST_DATA; i++) {
        target = &TEST_DATA[i];
        dim = target->dim;
        pattern = target->pattern;
        memcpy(matrix, target->matrix, sizeof(double)*dim*dim);

        // LU分解
        printf("=== dim: %d ===\n", dim);
        printf("Decomposition\n");
        ret = LU_decomposition(dim, pivot, matrix);
        if ((int32_t)LU_OK != ret) {
            goto EXIT_MAIN;
        }
        printf("Solve\n");
        // LU分解した結果を用いて演算
        for (j = 0; j < pattern; j++) {
            memcpy(exact, &(target->exacts[j * dim]), sizeof(double)*pattern);
            ret = mul_matvec(dim, target->matrix, exact, vec);
            if ((int32_t)RET_OK_MAIN != ret) {
                goto EXIT_MAIN;
            }
            ret = LU_solver(dim, pivot, matrix, vec);
            if ((int32_t)LU_OK != ret) {
                goto EXIT_MAIN;
            }
            ret = print_vec(dim, j, vec);
            if ((int32_t)RET_OK_MAIN != ret) {
                goto EXIT_MAIN;
            }
        }
        printf("\n");
    }
EXIT_MAIN:

    return 0;
}

static int32_t mul_matvec(int32_t dim, double *matrix, double *vec, double *out) {
    int32_t ret = (int32_t)RET_NG_MAIN;
    int row, col;
    double sum;

    if ((NULL != matrix) && (NULL != vec) && (NULL != out)) {
        for (row = 0; row < dim; row++) {
            sum = 0;

            for (col = 0; col < dim; col++) {
                sum += matrix[row * dim + col] * vec[col];
            }
            out[row] = sum;
        }
        ret = (int32_t)RET_OK_MAIN;
    }

    return ret;
}

static int32_t print_vec(int32_t dim, int32_t idx, double *out) {
    int32_t ret = (int32_t)RET_NG_MAIN;
    int32_t i;

    if (NULL != out) {
        printf(" [%d]", idx);
        for (i = 0; i < dim; i++) {
            printf(" %+.5f", out[i]);
        }
        printf("\n");
        ret = (int32_t)RET_OK_MAIN;
    }

    return ret;
}