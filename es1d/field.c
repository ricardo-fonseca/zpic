/*
 *  field.c
 *  zpic
 *
 *  Created by Ricardo Fonseca on 10/8/10.
 *  Copyright 2010 Centro de Física dos Plasmas. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>

#include "field.h"
#include "zdf.h"
#include "timer.h"

static double _field_time = 0.0;

double field_time()
{
	return _field_time;
}

/*********************************************************************************************
 
 Constructor / Destructor
 
 *********************************************************************************************/


void field_new( t_field *field, int nx, t_fld box, const float dt, 
	t_fftr_cfg *fft_backward )
{   	

	if ( nx % 2 ) {
		fprintf(stderr,"(*error*) Only even grid sizes are supported.");
		exit(-1);
	}

	// Number of guard cells for linear interpolation
	unsigned int gc[2] = {0,1}; 

	// Store pointer to required FFT configuration
	field -> fft_backward = fft_backward;

	// Allocate arrays
	scalar_grid_init( &field->E, nx, gc );
	cscalar_grid_init( &field->fE, nx/2+1, NULL );
	
	// zero fields
	scalar_grid_zero( &field->E );
	cscalar_grid_zero( &field->fE );
		
	// Set cell sizes and box limits
	field -> box = box;
	field -> dx = box / nx;

	// Set time step
	field -> dt = dt;

	// Reset iteration number
	field -> iter = 0;
	
}

void field_delete( t_field *field )
{
	scalar_grid_cleanup( &field->E );
	cscalar_grid_cleanup( &field->fE );
}

/*********************************************************************************************
 
 Diagnostics
 
 *********************************************************************************************/


void field_report( const t_field *field )
{
	char fname[] = "E1";
	
	float *buf = field -> E.s;

    t_zdf_grid_axis axis[1];
    axis[0] = (t_zdf_grid_axis) {
    	.min = 0.0,
    	.max = field->box,
    	.label = "x_1",
    	.units = "c/\\omega_p"
    };

    t_zdf_grid_info info = {
    	.ndims = 1,
    	.label = fname,
    	.units = "m_e c \\omega_p e^{-1}",
    	.axis = axis
    };

    info.nx[0] = field->E.nx;

    t_zdf_iteration iter = {
    	.n = field->iter,
    	.t = field -> iter * field -> dt,
    	.time_units = "1/\\omega_p"
    };

	zdf_save_grid( buf, &info, &iter, "field" );	
		
}


/*********************************************************************************************
 
 Field solver
 
 *********************************************************************************************/

void update_fE( t_field *field, const t_charge *charge )
{
    float complex * const restrict frho = charge -> frho.s;

    float complex * const restrict fEx = field -> fE.s;

    const float dk = fft_dk( field->E.nx, field->dx );
	unsigned int i;

	fEx[0] = 0;

    for( i = 1; i < field -> fE.nx; i++) {
    	float kx = i * dk;
    	fEx[i] = - I * frho[i] / kx;
    }

}


// This code operates with periodic boundaries
void field_update( t_field *field )
{
	int i;

	// Update E field
   	fftr_c2r( field -> fft_backward, field -> fE.s, field ->E.s );

   	// Update guard cells
    float* const restrict Ex = field -> E.s;
	
    const unsigned nx = field -> E.nx;

	// lower
	for (i = - field->E.gc[0]; i<0; i++) {
		Ex[ i ] = Ex[ nx + i ];
	}

	// upper
	for (i=0; i<field->E.gc[1]; i++) {
		Ex[ nx + i ] = Ex[ i ];
	}
	
}

void field_advance( t_field *field, const t_charge *charge )
{
	uint64_t t0 = timer_ticks();
	
	// Calculate fE
	update_fE( field, charge );

	// Update (real) E (also updates guard cells)
	field_update( field );

	// Advance internal iteration number
    field -> iter += 1;
	
    // Update timing information
	_field_time += timer_interval_seconds(t0, timer_ticks());
}
