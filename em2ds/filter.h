#ifndef __FILTER__
#define __FILTER__

#include "grid2d.h"

/**
 * @brief 2D low pass filter of scalar r2c data
 * 
 * Input data is expected to be the output of `fftr2d_r2c()`. Filtering is
 * done inplace.
 * 
 * @param grid 		scalar r2c data to filter
 * @param cutoff 	Cutoff frequency [kcx, kcy] as a fraction of the
 * 					Nyquist frequency
 */
void cscalar2d_r2c_filter( t_cscalar_grid2d * grid, const float cutoff[] );

/**
 * @brief 2D low pass filter of float3 r2c data
 * 
 * @param grid 		scalar r2c data to filter
 * @param cutoff	Cutoff frequency [kcx, kcy] as a fraction of the
 * 					Nyquist frequency
 */
void cfloat32d_r2c_filter( t_cfloat3_grid2d * grid, const float cutoff[] );

#endif
