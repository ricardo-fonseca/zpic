#include <stdlib.h>
#include <stdio.h>

#include "charge.h"
#include "zdf.h"
#include "fft.h"


void charge_new( t_charge *charge, int nx, t_fld box, float dt, t_fftr_cfg *fft )
{   

	if ( nx % 2 ) {
		fprintf(stderr,"(*error*) Only even grid sizes are supported.");
		exit(-1);
	}

	// Number of guard cells for linear interpolation
	unsigned int gc[2] = {1,2}; 
	
	// Store pointer to required FFT configuration
	charge -> fft_forward_2 = &fft[FORWARD_2];

	// Initialize grids
	scalar_grid_init( &charge->rho,   nx, gc );
	cscalar_grid_init( &charge->frho, nx+1, NULL );
	
	// Set cell sizes and box limits
	charge -> box = box;
	charge -> dx  = box / nx;

	// Initialize time information
	charge -> iter = 0;
	charge -> dt = dt;

    // Zero initial charge
    // This is only relevant for diagnostics, charge is always zeroed before deposition
	scalar_grid_zero( &charge -> rho );
	
	// Default is not to have a neutralizing background
	charge -> neutral.nx = 0;
}

void charge_init_neutral_bkg( t_charge *charge )
{
	scalar_grid_init( &charge->neutral, charge ->rho.nx,
					  (unsigned int *) charge ->rho.gc );
	scalar_grid_zero( &charge -> neutral );
}

void charge_delete( t_charge *charge )
{
	scalar_grid_cleanup( &charge -> rho );
	cscalar_grid_cleanup( &charge -> frho );

	if ( charge -> neutral.nx > 0 ) {
		scalar_grid_cleanup( &charge -> neutral );
	}
	
}

void charge_zero( t_charge *charge )
{
	scalar_grid_zero( &charge -> rho );
}

void charge_update( t_charge *charge )
{
	int i;

	float* restrict const rho = charge -> rho.s;
	float complex * restrict const frho = charge -> frho.s;
	
	int const gc0 = charge -> rho.gc[0];
	int const gc1 = charge -> rho.gc[1];
	int const nx  = charge -> rho.nx;

	// x
		
	// lower - add the values from upper boundary ( both gc and inside box )
	for (i=-gc0; i<gc1; i++) {
		rho[ i ] += rho[ nx + i ];
	}
	
	// upper - just copy the values from the lower boundary 
	for (i=-gc0; i<gc1; i++) {
		rho[ nx + i ] = rho[ i ];
	}

	// Add neutralizing background
	// This is preferable to initializing rho to this value before charge deposition
	// because it leads to less roundoff errors
	if ( charge -> neutral.nx > 0 ) {
		float* restrict const neutral = charge -> neutral.s;
		for(i=-gc0; i<nx+gc1; i++) {
			rho[i] += neutral[i];
		}
	}

	// Calculate frho
	// fftr_r2c( charge->fft_forward, rho, frho );

	// Calculate frho using DCT-II equivalent
	// even boundaries at -1/2 and nx-1/2

	float * rho_bnd = malloc( 2 * nx * sizeof(float) );

	for( i = 0; i < nx; i++ ) {
		rho_bnd[i] = rho[i];
	}

	for( i = 0; i < nx; i++ ) {
		rho_bnd[2*nx - 1 - i] = +rho[i];
	}

	fftr_r2c( charge->fft_forward_2, rho_bnd, frho );
	free(rho_bnd);

	// Filter charge
	for ( i = charge -> frho.nx/2; i < charge -> frho.nx; i++) {
		charge -> frho.s[i] = 0;
	}

	charge -> iter++;
}

void charge_update_neutral_bkg( t_charge *charge )
{
	int i;

	int const nx = charge->neutral.nx;
	int const gc0 = charge->neutral.gc[0];
	int const gc1 = charge->neutral.gc[1];

	// Update boundaries
	for (i = -gc0; i < gc1; i++) {
		charge->neutral.s[ i ] += charge->neutral.s[ nx + i ];
	}
	for (i = -gc0; i < gc1; i++) {
		charge->neutral.s[ nx + i ] = charge->neutral.s[ i ];
	}

	// Change sign
	for (i = -gc0; i < nx + gc1; i++) {
		charge->neutral.s[ i ] = -charge->neutral.s[ i ];
	}


}


void charge_report( const t_charge *charge )
{
	char vfname[] = "charge density";
		
	float *buf = charge -> rho.s;
	
    t_zdf_grid_axis axis[1];
    axis[0] = (t_zdf_grid_axis) {
    	.min = 0.0,
    	.max = charge->box,
    	.label = "x_1",
    	.units = "c/\\omega_p"
    };

    t_zdf_grid_info info = {
    	.ndims = 1,
    	.label = vfname,
    	.units = "e \\omega_p^2 / c",
    	.axis = axis
    };

    info.nx[0] = charge->rho.nx;

    t_zdf_iteration iter = {
    	.n = charge->iter,
    	.t = charge -> iter * charge -> dt,
    	.time_units = "1/\\omega_p"
    };

	zdf_save_grid( buf, &info, &iter, "CHARGE" );
		
}
