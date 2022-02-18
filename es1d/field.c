/**
 * @file field.c
 * @author Ricardo Fonseca
 * @brief Electric Field
 * @version 0.2
 * @date 2022-02-18
 * 
 * @copyright Copyright (c) 2022
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

/// Time spent advancing the electric field
static double _field_time = 0.0;

/**
 * @brief Time spent advancing the electric fields
 * 
 * @return      Time spent in seconds
 */
double field_time()
{
	return _field_time;
}

/*********************************************************************************************
 
 Constructor / Destructor
 
 *********************************************************************************************/

/**
 * @brief Initializes the electric field object
 * 
 * @param field 			Electric field object
 * @param nx 				Number of cells
 * @param box 				Physical box size
 * @param dt 				Simulation time step
 * @param fft_backward 		FFT configuration (shared with other objects)
 */
void field_new( t_field *field, int nx, float box, const float dt, 
	t_fftr_cfg *fft_backward )
{
	if ( nx % 2 ) {
		fprintf(stderr,"(*error*) Only even grid sizes are supported.");
		exit(-1);
	}

	// Number of guard cells for linear interpolation
	int gc[2] = {0,1}; 

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

/**
 * @brief Frees dynamic memory from electric field
 * 
 * @param field 	Electric field
 */
void field_delete( t_field *field )
{
	scalar_grid_cleanup( &field->E );
	cscalar_grid_cleanup( &field->fE );
}

/*********************************************************************************************
 
 Diagnostics
 
 *********************************************************************************************/

/**
 * @brief Saves electric field diagnostic information to disk
 * 
 * Field will be save in directory "field", guard cell values are
 * discarded.
 * 
 * @param field 	Electric field
 */
void field_report( const t_field *field )
{
	char fname[] = "E1";
	char flabel[] = "E_x";
	
	float *buf = field -> E.s;

    t_zdf_grid_axis axis[1];
    axis[0] = (t_zdf_grid_axis) {
    	.min = 0.0,
    	.max = field->box,
		.name = "x",
    	.label = "x",
    	.units = "c/\\omega_p"
    };

    t_zdf_grid_info info = {
    	.ndims = 1,
		.name = fname,
    	.label = flabel,
    	.units = "m_e c \\omega_p e^{-1}",
    	.axis = axis
    };

    info.count[0] = field->E.nx;

    t_zdf_iteration iter = {
		.name = "ITERATION",
    	.n = field->iter,
    	.t = field -> iter * field -> dt,
    	.time_units = "1/\\omega_p"
    };

	zdf_save_grid( (float *) buf, zdf_float32, &info, &iter, "field" );	
}


/*********************************************************************************************
 
 Field solver
 
 *********************************************************************************************/

/**
 * @brief Updates electric field from charge density
 * 
 * Note: We chose to set the k=0 component to 0
 * 
 * @param field 	Electric field
 * @param charge 	Charge density
 */
void update_fE( t_field *field, const t_charge *charge )
{
    float complex * const restrict frho = charge -> frho.s;

    float complex * const restrict fEx = field -> fE.s;

    const float dk = fft_dk( field->E.nx, field->dx );

	fEx[0] = 0;

    for( int i = 1; i < field -> fE.nx; i++) {
    	float kx = i * dk;
    	fEx[i] = - I * frho[i] / kx;
    }

}


/**
 * @brief Updates electric field from Fourier transform
 * 
 * It will also update guard cell values for field interpolation.
 * 
 * @param field 
 */
void field_update( t_field *field )
{
	// Update E field
   	fftr_c2r( field -> fft_backward, field -> fE.s, field ->E.s );

   	// Update guard cells
    float* const restrict Ex = field -> E.s;
	
    const unsigned nx = field -> E.nx;

	// lower
	for (int i = - field->E.gc[0]; i<0; i++) {
		Ex[ i ] = Ex[ nx + i ];
	}

	// upper
	for (int i=0; i<field->E.gc[1]; i++) {
		Ex[ nx + i ] = Ex[ i ];
	}
	
}

/**
 * @brief Advance electric field 1 timestep
 * 
 * Field is updated from the charge density. The routine will also update
 * guard cell values.
 * 
 * @param field 
 * @param charge 
 */
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
