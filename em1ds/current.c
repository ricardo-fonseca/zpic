/**
 * @file current.c
 * @author Ricardo Fonseca
 * @brief Electric current density
 * @version 0.2
 * @date 2022-02-21
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include "current.h"

#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>

#include "zdf.h"

/**
 * @brief Initializes Electric current density object
 * 
 * @param current 		Electric current density
 * @param nx 			Number of cells
 * @param box 			Physical box size
 * @param dt 			Simulation time step
 * @param fft_forward 	FFT configuration for transforming J to fJ
 * 						(shared with other objects)
 * @param filter 		Spectral filtering parameters
 */
void current_new( t_current *current, int nx, float box, float dt, 
	t_fftr_cfg *fft_forward, t_filter *filter )
{
	if ( nx % 2 ) {
		fprintf(stderr,"(*error*) Only even grid sizes are supported.");
		exit(-1);
	}

	// Number of guard cells for linear interpolation
	int gc[2] = {1,2};

	// Store pointer to required FFT configuration
	current -> fft_forward = fft_forward;

	// Store pointer to spectral filter data
	current -> filter = filter;

	// Initialize grids
	float3_grid_init( &current->J, nx, gc );
	cfloat3_grid_init( &current->fJ, nx/2+1, NULL );

	// Set cell sizes and box limits
	current -> box = box;
	current -> dx  = box / nx;

	// Initialize time information
	current -> iter = 0;
	current -> dt = dt;

    // Zero initial current

    // This is only relevant for diagnostics, current is always zeroed before deposition
	float3_grid_zero( &current -> J );

	// This avoids setting the fJ.x to zero in current update
	cfloat3_grid_zero( &current -> fJ );
}

/**
 * @brief Frees dynamic memory from electric current density
 * 
 * @param current Electric current density object
 */
void current_delete( t_current *current )
{
	float3_grid_cleanup( &current -> J );
	cfloat3_grid_cleanup( &current -> fJ );
}

/**
 * @brief Sets all electric current density values to zero
 * 
 * @param current Electric current density object
 */
void current_zero( t_current *current )
{
	float3_grid_zero( &current -> J );
}

/**
 * @brief Advances electric current density 1 time step
 * 
 * The routine will:
 * 1. Update the guard cells
 * 2. Calculate the transverse component of J (not needed in 1D)
 * 3. Get the Fourier transform of the current
 * 4. Apply spectral filtering (if configured)
 * 
 * @param current Electric current density object
 */
void current_update( t_current *current )
{
	float* restrict const Jx = current -> J.x;
	float* restrict const Jy = current -> J.y;
	float* restrict const Jz = current -> J.z;

	int const gc0 = current -> J.gc[0];
	int const gc1 = current -> J.gc[1];
	int const nx  = current -> J.nx;

	// x

	// lower - add the values from upper boundary ( both gc and inside box )
	for (int i=-gc0; i<gc1; i++) {
		Jx[ i ] += Jx[ nx + i ];
		Jy[ i ] += Jy[ nx + i ];
		Jz[ i ] += Jz[ nx + i ];
	}

	// upper - just copy the values from the lower boundary
	for (int i=-gc0; i<gc1; i++) {
		Jx[ nx + i ] = Jx[ i ];
		Jy[ nx + i ] = Jy[ i ];
		Jz[ nx + i ] = Jz[ i ];
	}

	// Calculate fJt

	// In 1D we have
	// fJt.x = 0; fJt.y = fJ.y; fJt.z = fJ.z

	// fJt.x will be always 0 (set on initialization) so there is no need to do an FFT

	// fftr_r2c( &current -> fft_forward, current -> J.x, current -> fJt.x );
	fftr_r2c( current -> fft_forward, Jy, current -> fJ.y );
	fftr_r2c( current -> fft_forward, Jz, current -> fJ.z );

	// Filter current
	if ( current -> filter -> type > FILTER_NONE ) {
	    float * const restrict Sk = current -> filter -> Sk;
	    complex float * const restrict fJtx = current -> fJ.x;
	    complex float * const restrict fJty = current -> fJ.y;
	    complex float * const restrict fJtz = current -> fJ.z;

		for ( int i = 0; i < current -> fJ.nx; i++) {
			fJtx[i] *= Sk[i];
			fJty[i] *= Sk[i];
			fJtz[i] *= Sk[i];
		}
	}

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
void current_report( const t_current *current, const int jc )
{
	float *buf = NULL;
	
	switch (jc) {
		case 0:
			buf = current -> J.x;
			break;
		case 1:
			buf = current -> J.y;
			break;
		case 2:
			buf = current -> J.y;
			break;
	}

	char vfname[16];	// Dataset name
	char vflabel[16];	// Dataset label (for plots)

    snprintf( vfname, 3, "J%1d", jc );
	char comp[] = {'x','y','z'};
    snprintf(vflabel,4,"J_%c",comp[jc]);
	
    t_zdf_grid_axis axis[1];
    axis[0] = (t_zdf_grid_axis) {
    	.min = 0.0,
    	.max = current->box,
		.name = "x",
    	.label = "x",
    	.units = "c/\\omega_p"
    };

    t_zdf_grid_info info = {
    	.ndims = 1,
        .name = vfname,
        .label = vflabel,
    	.units = "e \\omega_p^2 / c",
    	.axis = axis
    };

    info.count[0] = current->J.nx;

    t_zdf_iteration iter = {
        .name = "ITERATION",
    	.n = current->iter,
    	.t = current -> iter * current -> dt,
    	.time_units = "1/\\omega_p"
    };

    zdf_save_grid( (void *) buf, zdf_float32, &info, &iter, "CURRENT" );
}

