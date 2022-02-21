/**
 * @file filter.c
 * @author Ricardo Fonseca
 * @brief 2D low pass spectral filtering of real data
 * @version 0.1
 * @date 2022-02-17
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#include "grid2d.h"
#include <stdlib.h>

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
void cscalar2d_r2c_filter( t_cscalar_grid2d * grid, const float cutoff[] )
{
	const int nkx = grid -> nx[1];
	const int nky = grid -> nx[0];
	const int nrow = grid -> nrow;

	// Cutoff value in kx
	const int kcx = cutoff[0] * ( nkx - 1 );

	// Cutoff value in ky
	const int kcy = cutoff[0] * ( grid -> nx[0] / 2 );

	for ( int i = 0; i < nkx; i++) {
		const int kx = i;
		for( int j = 0; j < nky; j++) {
			const int ky = abs( (j <= nky/2) ? j : ( j - nky ) );
			if (( kx > kcx ) || (ky > kcy)) {
				grid -> s[ i*nrow + j ] = 0;
			}
		}
	}

}

/**
 * @brief 2D low pass filter of float3 r2c data
 * 
 * @param grid 		scalar r2c data to filter
 * @param cutoff	Cutoff frequency [kcx, kcy] as a fraction of the
 * 					Nyquist frequency
 */
void cfloat32d_r2c_filter( t_cfloat3_grid2d * grid, const float cutoff[] )
{
	const int nkx = grid -> nx[1];
	const int nky = grid -> nx[0];
	const int nrow = grid -> nrow;

	// Cutoff value in kx
	const int kcx = cutoff[0] * ( nkx - 1 );

	// Cutoff value in ky
	const int kcy = cutoff[0] * ( grid -> nx[0] / 2 );

	for ( int i = 0; i < nkx; i++) {
		const int kx = i;
		for( int j = 0; j < nky; j++) {
			const int ky = abs( (j <= nky/2) ? j : ( j - nky ) );
			if (( kx > kcx ) || (ky > kcy)) {
				grid -> x[ i*nrow + j ] = 0;
				grid -> y[ i*nrow + j ] = 0;
				grid -> z[ i*nrow + j ] = 0;
			}
		}
	}

}
