#include "grid2d.h"
#include <stdlib.h>

void cscalar2d_r2c_filter( t_cscalar_grid2d * grid, const float cutoff[] )
{
	int i, j;

	const int nkx = grid -> nx[1];
	const int nky = grid -> nx[0];
	const int nrow = grid -> nrow;

	// Cutoff value in kx
	const int kcx = cutoff[0] * ( nkx - 1 );

	// Cutoff value in ky
	const int kcy = cutoff[0] * ( grid -> nx[0] / 2 );

	for ( i = 0; i < nkx; i++) {
		const int kx = i;
		for( j = 0; j < nky; j++) {
			const int ky = abs( (j <= nky/2) ? j : ( j - nky ) );
			if (( kx > kcx ) || (ky > kcy)) {
				grid -> s[ i*nrow + j ] = 0;
			}
		}
	}

}



void cvfld2d_r2c_filter( t_cvfld_grid2d * grid, const float cutoff[] )
{
	int i, j;

	const int nkx = grid -> nx[1];
	const int nky = grid -> nx[0];
	const int nrow = grid -> nrow;

	// Cutoff value in kx
	const int kcx = cutoff[0] * ( nkx - 1 );

	// Cutoff value in ky
	const int kcy = cutoff[0] * ( grid -> nx[0] / 2 );

	for ( i = 0; i < nkx; i++) {
		const int kx = i;
		for( j = 0; j < nky; j++) {
			const int ky = abs( (j <= nky/2) ? j : ( j - nky ) );
			if (( kx > kcx ) || (ky > kcy)) {
				grid -> x[ i*nrow + j ] = 0;
				grid -> y[ i*nrow + j ] = 0;
				grid -> z[ i*nrow + j ] = 0;
			}
		}
	}

}
