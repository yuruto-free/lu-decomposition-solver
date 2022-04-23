#ifndef LU_DECOMP_SOLVER_H__
#define LU_DECOMP_SOLVER_H__

#include <stdint.h>

#define LU_OK (0)
#define LU_NG (1)

/**
 * @brief LU分解
 * 
 * @param[in]    dim    次元数
 * @param[inout] pivot  ピボット
 *                      サイズ：dim
 * @param[inout] matrix 係数行列
 *                      サイズ：dim * dim
 * @return LU_OK  正常終了
 * @return LU_NG  異常終了
 */
int32_t LU_decomposition(int32_t dim, int32_t *pivot, double *matrix);

/**
 * @brief LU分解を用いた連立一次方程式の求解
 * @param[in]    dim    次元数
 * @param[in]    pivot  ピボット
 *                      サイズ：dim
 * @param[in]    matrix LU分解済みの係数行列
 *                      サイズ：dim * dim
 * @param[inout] vec    右辺ベクトル
 *                      サイズ：dim
 * @return LU_OK  正常終了
 * @return LU_NG  異常終了
 */
int32_t LU_solver(int32_t dim, int32_t *pivot, double *matrix, double *vec);

#endif