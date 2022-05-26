/**
 * @file current.c
 * @author Ricardo Fonseca
 * @brief Electric current density
 * @version 0.2
 * @date 2022-02-04
 * 
 * @copyright Copyright (c) 2022
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

/**
 * @brief Initializes Electric current density object
 * 
 * @param current 		Electric current density
 * @param nx 			Number of cells
 * @param box 			Physical box size
 * @param dt 			Simulation time step
 */
void current_new( t_current *current, const int nx[], float box[], float dt )
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
	float3_grid2d_init( &current->J, (unsigned int *) nx, gc );
	cfloat3_grid2d_init( &current->fJ, fnx, NULL );

	// Initializ FFT transform
	fftr2d_init_cfg( &current -> fft_forward, nx[0], nx[1],
	    current -> J.nrow, FFT_FORWARD );

	fftr2d_init_cfg( &current -> fft_backward, nx[0], nx[1],
	    current -> J.nrow, FFT_BACKWARD );

	// Set cell sizes and box limits
	for(int i = 0; i<2; i++){
		current -> box[i] = box[i];
		current -> dx[i]  = box[i] / nx[i];
	}

	// Initialize time information
	current -> iter = 0;
	current -> dt = dt;

  // Zero initial current

  // This is only relevant for diagnostics, current is always zeroed before deposition
	float3_grid2d_zero( &current -> J );

}

/**
 * @brief Frees dynamic memory from electric current density
 * 
 * @param current Electric current density object
 */
void current_delete( t_current *current )
{
	float3_grid2d_cleanup( &current -> J );
	cfloat3_grid2d_cleanup( &current -> fJ );

	fftr2d_cleanup_cfg( &current -> fft_forward );

}

/**
 * @brief Sets all electric current density values to zero
 * 
 * @param current Electric current density object
 */
void current_zero( t_current *current )
{
	float3_grid2d_zero( &current -> J );
}

/**
 * @brief Advances electric current density 1 time step
 * 
 * The routine will:
 * 1. Update the guard cells
 * 2. Calculate the transverse component of J (not needed in 1D)
 * 3. Get the Fourier transform of the current
 * 4. Apply spectral filtering
 * 
 * @param current Electric current density object
 */
void current_update( t_current *current )
{
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
	for (int j = -gc10; j < nx1 + gc11; j++) {

		// lower - add the values from upper boundary ( both gc and inside box )
		for (int i = -gc00; i < gc01; i++) {
			Jx[ i + j*nrow ] += Jx[ nx0 + i + j*nrow ];
			Jy[ i + j*nrow ] += Jy[ nx0 + i + j*nrow ];
			Jz[ i + j*nrow ] += Jz[ nx0 + i + j*nrow ];
		}

		// upper - just copy the values from the lower boundary
		for (int i = -gc00; i < gc01; i++) {
			Jx[ nx0 + i + j*nrow ] = Jx[ i + j*nrow ];
			Jy[ nx0 + i + j*nrow ] = Jy[ i + j*nrow ];
			Jz[ nx0 + i + j*nrow ] = Jz[ i + j*nrow ];
		}

	}

	// y
	for (int i = -gc00; i < nx0 + gc01; i++) {

		// lower - add the values from upper boundary ( both gc and inside box )
		for (int j=-gc10; j<gc11; j++) {
			Jx[ i + j*nrow ] += Jx[ i + (nx1+j)*nrow ];
			Jy[ i + j*nrow ] += Jy[ i + (nx1+j)*nrow ];
			Jz[ i + j*nrow ] += Jz[ i + (nx1+j)*nrow ];
		}

		// upper - just copy the values from the lower boundary
		for (int j=-gc10; j<gc11; j++) {
			Jx[ i + (nx1+j)*nrow ] = Jx[ i + j*nrow ];
			Jy[ i + (nx1+j)*nrow ] = Jy[ i + j*nrow ];
			Jz[ i + (nx1+j)*nrow ] = Jz[ i + j*nrow ];
		}

	}

	// Calculate fJ
	float complex* restrict const fJx = current -> fJ.x;
	float complex* restrict const fJy = current -> fJ.y;
	float complex* restrict const fJz = current -> fJ.z;

	fftr2d_r2c( &current -> fft_forward, Jx, fJx );
	fftr2d_r2c( &current -> fft_forward, Jy, fJy );
	fftr2d_r2c( &current -> fft_forward, Jz, fJz );

	// Filter current
	const float cutoff[2] = {0.5f,0.5f};
	cfloat32d_r2c_filter( &current -> fJ, cutoff );

	// Advance iteration counter
	current -> iter++;

}

/**
 * @brief Saves electric current density diagnostic information to disk
 * 
 * Saves the selected current density component to disk in directory
 * "CURRENT". Guard cell values are discarded.
 * 
 * @param current Electric current object
 * @param jc Current component to save, must be one of {0,1,2}
 */
void current_report( const t_current *current, const char jc )
{
	float *f, *buf, *p;
	char comp;

	switch (jc) {
		case 0:
			f = current->J.x;
			comp = 'x';
			break;
		case 1:
			f = current->J.y;
			comp = 'y';
			break;
		case 2:
			f = current->J.z;
			comp = 'z';
			break;
		default:
			fprintf(stderr,"(*error*) Invalid current component selected, returning\n");
			return;
	}

	char vfname[16];	// Dataset name
	char vflabel[16];	// Dataset label (for plots)

    snprintf( vfname, 3, "J%1u", jc );
    snprintf(vflabel,4,"J_%c",comp);

	// Pack the information
	buf = malloc( current->J.nx[0]*current->J.nx[1]*sizeof(float) );
    p = buf;

	for( int j = 0; j < current->J.nx[1]; j++) {
		for ( int i = 0; i < current->J.nx[0]; i++ ) {
			p[i] = f[i];
		}
		p += current->J.nx[0];
		f += current->J.nrow;
	}

  t_zdf_grid_axis axis[2];
  axis[0] = (t_zdf_grid_axis) {
  	.min = 0.0,
  	.max = current->box[0],
  	.name = "x",
	.label = "x",
  	.units = "c/\\omega_p"
  };

  axis[1] = (t_zdf_grid_axis) {
  	.min = 0.0,
  	.max = current->box[1],
  	.name = "y",
	.label = "y",
  	.units = "c/\\omega_p"
  };

  t_zdf_grid_info info = {
  	.ndims = 2,
  	.name = vfname,
	.label = vflabel,
  	.units = "e \\omega_p^2 / c",
  	.axis = axis
  };

  info.count[0] = current->J.nx[0];
  info.count[1] = current->J.nx[1];

  t_zdf_iteration iter = {
	.name = "ITERATION",
  	.n = current->iter,
  	.t = current -> iter * current -> dt,
  	.time_units = "1/\\omega_p"
  };

	zdf_save_grid( (void *) buf, zdf_float32, &info, &iter, "CURRENT" );

	// free local data
	free( buf );

}
