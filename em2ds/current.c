/*
 *  current.c
 *  zpic
 *
 *  Created by Ricardo Fonseca on 12/8/10.
 *  Copyright 2010 Centro de Física dos Plasmas. All rights reserved.
 *
 */

#include "current.h"

#include <math.h>

#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "fft.h"
#include "filter.h"
#include "zdf.h"

void current_new( t_current *current, const unsigned int nx[], t_fld box[], float dt )
{

	if ( nx[0] % 2 || nx[1] % 2 ) {
		fprintf(stderr,"(*error*) Only even grid sizes are supported.");
		exit(-1);
	}

	// Number of guard cells for linear interpolation
	const unsigned int gc[2][2] = {{1,2},
		            {1,2}};

	// The FFT result will be transposed
	const unsigned int fnx[2] = { nx[1], nx[0]/2+1 };

	// Initialize grids
	vfld_grid2d_init( &current->J, nx, gc );
	cvfld_grid2d_init( &current->fJt, fnx, NULL );

	// Initializ FFT transform
	fftr2d_init_cfg( &current -> fft_forward, nx[0], nx[1],
	    current -> J.nrow, FFT_FORWARD );

	fftr2d_init_cfg( &current -> fft_backward, nx[0], nx[1],
	    current -> J.nrow, FFT_BACKWARD );

	// Set cell sizes and box limits
	unsigned i;
	for(i = 0; i<2; i++){
		current -> box[i] = box[i];
		current -> dx[i]  = box[i] / nx[i];
	}

	// Initialize time information
	current -> iter = 0;
	current -> dt = dt;

  // Zero initial current

  // This is only relevant for diagnostics, current is always zeroed before deposition
	vfld_grid2d_zero( &current -> J );

	cvfld_grid2d_zero( &current -> fJt );

}

void current_delete( t_current *current )
{
	vfld_grid2d_cleanup( &current -> J );
	cvfld_grid2d_cleanup( &current -> fJt );

	fftr2d_cleanup_cfg( &current -> fft_forward );

}

void current_zero( t_current *current )
{
	vfld_grid2d_zero( &current -> J );
}

void current_update( t_current *current )
{
	int i,j;
	const int nrow = current->J.nrow;
	float* restrict const Jx = current -> J.x;
	float* restrict const Jy = current -> J.y;
	float* restrict const Jz = current -> J.z;

	int const nx0 = current -> J.nx[0];
	int const nx1 = current -> J.nx[1];

	int const gc00 = current -> J.gc[0][0];
	int const gc01 = current -> J.gc[0][1];
	int const gc10 = current -> J.gc[1][0];
	int const gc11 = current -> J.gc[1][1];


	// x
	for (j = -gc10; j < nx1 + gc11; j++) {

		// lower - add the values from upper boundary ( both gc and inside box )
		for (i = -gc00; i < gc01; i++) {
			Jx[ i + j*nrow ] += Jx[ nx0 + i + j*nrow ];
			Jy[ i + j*nrow ] += Jy[ nx0 + i + j*nrow ];
			Jz[ i + j*nrow ] += Jz[ nx0 + i + j*nrow ];
		}

		// upper - just copy the values from the lower boundary
		for (i = -gc00; i < gc01; i++) {
			Jx[ nx0 + i + j*nrow ] = Jx[ i + j*nrow ];
			Jy[ nx0 + i + j*nrow ] = Jy[ i + j*nrow ];
			Jz[ nx0 + i + j*nrow ] = Jz[ i + j*nrow ];
		}

	}

	// y
	for (i = -gc00; i < nx0 + gc01; i++) {

		// lower - add the values from upper boundary ( both gc and inside box )
		for (j=-gc10; j<gc11; j++) {
			Jx[ i + j*nrow ] += Jx[ i + (nx1+j)*nrow ];
			Jy[ i + j*nrow ] += Jy[ i + (nx1+j)*nrow ];
			Jz[ i + j*nrow ] += Jz[ i + (nx1+j)*nrow ];
		}

		// upper - just copy the values from the lower boundary
		for (j=-gc10; j<gc11; j++) {
			Jx[ i + (nx1+j)*nrow ] = Jx[ i + j*nrow ];
			Jy[ i + (nx1+j)*nrow ] = Jy[ i + j*nrow ];
			Jz[ i + (nx1+j)*nrow ] = Jz[ i + j*nrow ];
		}

	}

	// Calculate fJt

	float complex* restrict const fJtx = current -> fJt.x;
	float complex* restrict const fJty = current -> fJt.y;
	float complex* restrict const fJtz = current -> fJt.z;
	const int fnrow = current->fJt.nrow;

	fftr2d_r2c( &current -> fft_forward, Jx, fJtx );
	fftr2d_r2c( &current -> fft_forward, Jy, fJty );
	fftr2d_r2c( &current -> fft_forward, Jz, fJtz );

	const float dkx = fft_dk( current->J.nx[0], current->dx[0] );
	const float dky = fft_dk( current->J.nx[1], current->dx[1] );

	// Calculates the transverse component of the current density
	// $\mathbf{j_T} = \mathbf{j} - \frac{\mathbf{k}( \mathbf{k} \cdot \mathbf{j})}{k^2}$

	fJtx[0] = 0;
	fJty[0] = 0;
	//fJtz[0] = 0;

	// i = 0
	for( j = 1; j < current -> fJt.nx[0]; j++) {
		float ky = ((j <= current -> fJt.nx[0]/2) ? j : ( j - (int) current -> fJt.nx[0] ) ) * dky;
		fJty[j] -= fJty[j]/ky;
	}


	for ( i = 1; i < current -> fJt.nx[1]; i++) {
		float kx = i * dkx;
		for( j = 0; j < current -> fJt.nx[0]; j++) {
			unsigned int idx =  i * fnrow + j ;

			float ky = ((j <= current -> fJt.nx[0]/2) ? j : ( j - (int) current -> fJt.nx[0] ) ) * dky;
			float ksqr = kx*kx + ky*ky;
			float complex kdJ = kx * fJtx[idx] + ky * fJty[idx];

			fJtx[idx] -= (kx * kdJ)/ksqr;
			fJty[idx] -= (ky * kdJ)/ksqr;
		}
	}

	// Filter current
	const float cutoff[2] = {0.5f,0.5f};
	cvfld2d_r2c_filter( &current -> fJt, cutoff );

	current -> iter++;

}

void current_report( const t_current *current, const char jc )
{
	float *f;
	float *buf, *p;
	unsigned int i, j;
	char vfname[3];

	// Pack the information
	buf = malloc( current->J.nx[0]*current->J.nx[1]*sizeof(float) );
    p = buf;

	vfname[0] = 'J';

	switch (jc) {
		case 0:
			f = current->J.x;
			vfname[1] = '1';
			break;
		case 1:
			f = current->J.y;
			vfname[1] = '2';
			break;
		case 2:
			f = current->J.z;
			vfname[1] = '3';
			break;
	}
	vfname[2] = 0;

	for( j = 0; j < current->J.nx[1]; j++) {
		for ( i = 0; i < current->J.nx[0]; i++ ) {
			p[i] = f[i];
		}
		p += current->J.nx[0];
		f += current->J.nrow;
	}

  t_zdf_grid_axis axis[2];
  axis[0] = (t_zdf_grid_axis) {
  	.min = 0.0,
  	.max = current->box[0],
  	.label = "x_1",
  	.units = "c/\\omega_p"
  };

  axis[1] = (t_zdf_grid_axis) {
  	.min = 0.0,
  	.max = current->box[1],
  	.label = "x_2",
  	.units = "c/\\omega_p"
  };

  t_zdf_grid_info info = {
  	.ndims = 2,
  	.label = vfname,
  	.units = "e \\omega_p^2 / c",
  	.axis = axis
  };

  info.nx[0] = current->J.nx[0];
  info.nx[1] = current->J.nx[1];

  t_zdf_iteration iter = {
  	.n = current->iter,
  	.t = current -> iter * current -> dt,
  	.time_units = "1/\\omega_p"
  };

	zdf_save_grid( buf, &info, &iter, "CURRENT" );

	// free local data
	free( buf );

}
