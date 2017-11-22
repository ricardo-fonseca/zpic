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

#include "zdf.h"

void current_new( t_current *current, int nx, t_fld box, float dt )
{   

	// Number of guard cells for linear interpolation
	int gc[2] = {1,2}; 
	
	// Allocate global array
	size_t size;
	
	size = gc[0] + nx + gc[1];
	
	current->J_buf = malloc( size * sizeof( t_vfld ) );
	assert( current->J_buf );

	// store nx and gc values
	current->nx = nx;
	current->gc[0] = gc[0];
	current->gc[1] = gc[1];
	
	// Make J point to cell [0]
	current->J = current->J_buf + gc[0];
	
	// Set cell sizes and box limits
	current -> box = box;
	current -> dx  = box / nx;

	// Clear smoothing options
	current -> smooth = (t_smooth) {
		.xtype = none,
		.xlevel = 0
	};

	// Initialize time information
	current -> iter = 0;
	current -> dt = dt;

	// Default to periodic boundaries
	current -> bc_type = CURRENT_BC_PERIODIC;

    // Zero initial current
    // This is only relevant for diagnostics, current is always zeroed before deposition
	current_zero( current );
	
}

void current_delete( t_current *current )
{
	free( current->J_buf );
	
	current->J_buf = NULL;
	
}

void current_zero( t_current *current )
{
	// zero fields
	size_t size;
	
	size = (current->gc[0] + current->nx + current->gc[1]) * sizeof( t_vfld );
	memset( current->J_buf, 0, size );
	
}

void current_update_gc( t_current *current )
{
	if ( current -> bc_type == CURRENT_BC_PERIODIC ) {
		t_vfld* restrict const J = current -> J;
		const int nx = current -> nx;

		// lower - add the values from upper boundary ( both gc and inside box )
		for (int i=-current->gc[0]; i<current->gc[1]; i++) {
			J[ i ].x += J[ nx + i ].x;
			J[ i ].y += J[ nx + i ].y;
			J[ i ].z += J[ nx + i ].z;
		}
		
		// upper - just copy the values from the lower boundary 
		for (int i=-current->gc[0]; i<current->gc[1]; i++) {
			J[ nx + i ].x = J[ i ].x;
			J[ nx + i ].y = J[ i ].y;
			J[ nx + i ].z = J[ i ].z;
		}
	}
}

void current_update( t_current *current )
{
	
	// Boundary conditions / guard cells
	current_update_gc( current );

	// Smoothing
	current_smooth( current );

	// Advance iteration number
	current -> iter++;
	
}

void current_report( const t_current *current, const char jc )
{
	t_vfld *f;
	int i;
	char vfname[3];
		
	// Pack the information
	float buf[current->nx];

	f = current->J;
	vfname[0] = 'J';
	
	switch (jc) {
		case 0:
			for ( i = 0; i < current->nx; i++ ) {
				buf[i] = f[i].x;
			}
			vfname[1] = '1';
			break;
		case 1:
			for ( i = 0; i < current->nx; i++ ) {
				buf[i] = f[i].y;
			}
			vfname[1] = '2';
			break;
		case 2:
			for ( i = 0; i < current->nx; i++ ) {
				buf[i] = f[i].z;
			}
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

    info.nx[0] = current->nx;

    t_zdf_iteration iter = {
    	.n = current->iter,
    	.t = current -> iter * current -> dt,
    	.time_units = "1/\\omega_p"
    };

	zdf_save_grid( buf, &info, &iter, "CURRENT" );
		
}

/*
 * get_smooth_comp
 *  Gets the value of the compensator kernel for an n pass binomial kernel
 */

void get_smooth_comp( int n, t_fld* sa, t_fld* sb) {
	t_fld a,b,total;

	a = -1;
	b = (4.0 + 2.0*n)/n;
	total = 2*a + b;

	*sa = a / total;
	*sb = b / total;
}


void kernel_x( t_current* const current, const t_fld sa, const t_fld sb ){

	int i;
	t_vfld* restrict const J = current -> J;


	t_vfld fl = J[-1];
	t_vfld f0 = J[ 0];

	for( i = 0; i < current -> nx; i++) {

		t_vfld fu = J[i + 1];

		t_vfld fs;

		fs.x = sa * fl.x + sb * f0.x + sa * fu.x;
		fs.y = sa * fl.y + sb * f0.y + sa * fu.y;
		fs.z = sa * fl.z + sb * f0.z + sa * fu.z;

		J[i] = fs;

		fl = f0;
		f0 = fu;

	}

	// Update x boundaries for periodic boundaries
	if ( current -> bc_type == CURRENT_BC_PERIODIC ) {
		for(i = -current->gc[0]; i<0; i++) 
			J[ i ] = J[ current->nx + i ];

		for (i=0; i<current->gc[1]; i++)
			J[ current->nx + i ] = J[ i ];
	}


}


void current_smooth( t_current* const current ) {

	// filter kernel [sa, sb, sa]
	t_fld sa, sb;

    // x-direction filtering
    if ( current -> smooth.xtype != none ) {
    	// binomial filter
    	sa = 0.25; sb = 0.5;
    	for( int i = 0; i < current -> smooth.xlevel; i++) {
    		kernel_x( current, 0.25, 0.5 );
    	}

    	// Compensator
    	if ( current -> smooth.xtype == compensated ) {
    		get_smooth_comp( current -> smooth.xlevel, &sa, &sb );
    		kernel_x( current, sa, sb );
    	}
    }

}

