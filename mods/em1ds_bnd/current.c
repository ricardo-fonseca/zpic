/*
 *  current.c
 *  zpic
 *
 *  Created by Ricardo Fonseca on 12/8/10.
 *  Copyright 2010 Centro de Física dos Plasmas. All rights reserved.
 *
 */

#include "current.h"

#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>

#include "zdf.h"

void current_new( t_current *current, int nx, t_fld box, float dt, t_fftr_cfg *fft )
{   

	if ( nx % 2 ) {
		fprintf(stderr,"(*error*) Only even grid sizes are supported.");
		exit(-1);
	}

	// Number of guard cells for linear interpolation
	unsigned int gc[2] = {1,2}; 
	
	// Store pointer to required FFT configuration
	current -> fft_forward = &fft[FORWARD];

	// Initialize grids
	vfld_grid_init( &current->J, nx, gc );
	cvfld_grid_init( &current->fJt, nx/2+1, NULL );
	
	// Set cell sizes and box limits
	current -> box = box;
	current -> dx  = box / nx;

	// Initialize time information
	current -> iter = 0;
	current -> dt = dt;

    // Zero initial current

    // This is only relevant for diagnostics, current is always zeroed before deposition
	vfld_grid_zero( &current -> J );

	// This avoids setting the fJt.x to zero in current update
	cvfld_grid_zero( &current -> fJt );
	
}

void current_delete( t_current *current )
{
	vfld_grid_cleanup( &current -> J );
	cvfld_grid_cleanup( &current -> fJt );
	
}

void current_zero( t_current *current )
{
	vfld_grid_zero( &current -> J );
}

void current_update( t_current *current )
{
	int i;

	float* restrict const Jx = current -> J.x;
	float* restrict const Jy = current -> J.y;
	float* restrict const Jz = current -> J.z;
	
	int const gc0 = current -> J.gc[0];
	int const gc1 = current -> J.gc[1];
	int const nx  = current -> J.nx;

	// x
		
	// lower - add the values from upper boundary ( both gc and inside box )
	for (i=-gc0; i<gc1; i++) {
		Jx[ i ] += Jx[ nx + i ];
		Jy[ i ] += Jy[ nx + i ];
		Jz[ i ] += Jz[ nx + i ];
	}
	
	// upper - just copy the values from the lower boundary 
	for (i=-gc0; i<gc1; i++) {
		Jx[ nx + i ] = Jx[ i ];
		Jy[ nx + i ] = Jy[ i ];
		Jz[ nx + i ] = Jz[ i ];
	}

	// Calculate fJt
	
	// In 1D we have
	// fJt.x = 0; fJt.y = fJ.y; fJt.z = fJ.z

	// fJt.x will be always 0 (set on initialization) so there is no need to do an FFT
	
	// fftr_r2c( &current -> fft_forward, current -> J.x, current -> fJt.x );
	fftr_r2c( current -> fft_forward, current -> J.y, current -> fJt.y );
	fftr_r2c( current -> fft_forward, current -> J.z, current -> fJt.z );

	// Filter current
	for ( i = current -> fJt.nx/2; i < current -> fJt.nx; i++) {
		current -> fJt.y[i] = 0;
		current -> fJt.z[i] = 0;
	}

	current -> iter++;
	
}

void current_report( const t_current *current, const char jc )
{
	char vfname[3];
		
	float *buf;

	vfname[0] = 'J';
	
	switch (jc) {
		case 0:
			buf = current -> J.x;
			vfname[1] = '1';
			break;
		case 1:
			buf = current -> J.y;
			vfname[1] = '2';
			break;
		case 2:
			buf = current -> J.y;
			vfname[1] = '3';
			break;
	}
	vfname[2] = 0;
	
    t_zdf_grid_axis axis[1];
    axis[0] = (t_zdf_grid_axis) {
    	.min = 0.0,
    	.max = current->box,
    	.label = "x_1",
    	.units = "c/\\omega_p"
    };

    t_zdf_grid_info info = {
    	.ndims = 1,
    	.label = vfname,
    	.units = "e \\omega_p^2 / c",
    	.axis = axis
    };

    info.nx[0] = current->J.nx;

    t_zdf_iteration iter = {
    	.n = current->iter,
    	.t = current -> iter * current -> dt,
    	.time_units = "1/\\omega_p"
    };

	zdf_save_grid( buf, &info, &iter, "CURRENT" );
		
}

