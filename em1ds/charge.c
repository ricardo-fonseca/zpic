#include <stdlib.h>
#include <stdio.h>

#include "charge.h"
#include "zdf.h"
#include "fft.h"
#include "filter.h"

/**
 * @brief Initializes Electric charge density object
 * 
 * @param charge 		Electric charge density
 * @param nx 			Number of cells
 * @param box 			Physical box size
 * @param dt 			Simulation time step
 * @param fft_forward 	FFT configuration for transforming rho to frho
 * 						(shared with other objects)
 * @param filter 		Spectral filtering parameterso
 */
void charge_new( t_charge *charge, int nx, float box, float dt, t_fftr_cfg *fft_forward,
                 t_filter *filter )
{

	if ( nx % 2 ) {
		fprintf(stderr,"(*error*) Only even grid sizes are supported.");
		exit(-1);
	}

	// Number of guard cells for linear interpolation
	int gc[2] = {1,2};

	// Store pointer to required FFT configuration
	charge -> fft_forward = fft_forward;

	// Store pointer to spectral filter data
	charge -> filter = filter;

	// Initialize grids
	scalar_grid_init( &charge->rho, nx, gc );
	cscalar_grid_init( &charge->frho, nx/2+1, NULL );

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

/**
 * @brief Initializes neutralizing background structures
 * 
 * @param charge 	Electric charge density
 */
void charge_init_neutral_bkg( t_charge *charge )
{
	scalar_grid_init( &charge->neutral, charge ->rho.nx, charge ->rho.gc );
	scalar_grid_zero( &charge -> neutral );
}

/**
 * @brief Frees dynamic memory from electric charge density
 * 
 * @param charge 	Electric charge density
 */
void charge_delete( t_charge *charge )
{
	scalar_grid_cleanup( &charge -> rho );
	cscalar_grid_cleanup( &charge -> frho );

	if ( charge -> neutral.nx > 0 ) {
		scalar_grid_cleanup( &charge -> neutral );
	}

}

/**
 * @brief Sets all electric charge density values to zero
 * 
 * @param charge 	Electric charge density
 */
void charge_zero( t_charge *charge )
{
	scalar_grid_zero( &charge -> rho );
}

/**
 * @brief Advances electric charge density 1 time step
 * 
 * The routine will:
 * 1. Update the guard cells
 * 2. Add neutralizing background (if configured)
 * 3. Get the fourier transform of the charge
 * 4. Apply spectral filtering (if configured)
 * 
 * @param charge 	Electric charge density
 */
void charge_update( t_charge *charge )
{
	float* restrict const rho = charge -> rho.s;
	float complex * restrict const frho = charge -> frho.s;

	int const gc0 = charge -> rho.gc[0];
	int const gc1 = charge -> rho.gc[1];
	int const nx  = charge -> rho.nx;

	// x
	// lower - add the values from upper boundary ( both gc and inside box )
	for (int i=-gc0; i<gc1; i++) {
		rho[ i ] += rho[ nx + i ];
	}

	// upper - just copy the values from the lower boundary
	for (int i=-gc0; i<gc1; i++) {
		rho[ nx + i ] = rho[ i ];
	}

	// Add neutralizing background
	// This is preferable to initializing rho to this value before charge deposition
	// because it leads to less roundoff errors
	if ( charge -> neutral.nx > 0 ) {
		float* restrict const neutral = charge -> neutral.s;
		for(int i=-gc0; i<nx+gc1; i++) {
			rho[i] += neutral[i];
		}
	}

	// Calculate frho
	fftr_r2c( charge->fft_forward, rho, frho );

	// Filter charge
	if ( charge -> filter -> type > FILTER_NONE ) {
	    float * const restrict Sk = charge -> filter -> Sk;
		for (int i = 0; i < charge -> frho.nx; i++) {
			charge -> frho.s[i] *= Sk[i];
		}
	}

	charge -> iter++;
}

/**
 * @brief Updates neutralizing background values
 * 
 * The routine expects the charge -> neutral grid to have the charge
 * density that needs to be neutralized. The routine will update guard
 * cell values and reverse the sign of the charge.
 * 
 * @param charge 
 */
void charge_update_neutral_bkg( t_charge *charge )
{
	int const nx  = charge->neutral.nx;
	int const gc0 = charge->neutral.gc[0];
	int const gc1 = charge->neutral.gc[1];

	// Update boundaries
	for (int i = -gc0; i < gc1; i++) {
		charge->neutral.s[ i ] += charge->neutral.s[ nx + i ];
	}

	for (int i =  -gc0; i < gc1; i++) {
		charge->neutral.s[ nx + i ] = charge->neutral.s[ i ];
	}

	// Change sign
	for (int i = -gc0; i < nx + gc1; i++) {
		charge->neutral.s[ i ] = -charge->neutral.s[ i ];
	}

}


void charge_report( const t_charge *charge )
{
	char vfname[] = "charge_density";
	char vflabel[] = "\\rho";

	float *buf = charge -> rho.s;

    t_zdf_grid_axis axis[1];
    axis[0] = (t_zdf_grid_axis) {
    	.min = 0.0,
    	.max = charge->box,
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

    info.count[0] = charge->rho.nx;

    t_zdf_iteration iter = {
		.name = "ITERATION",
    	.n = charge->iter,
    	.t = charge -> iter * charge -> dt,
    	.time_units = "1/\\omega_p"
    };

	zdf_save_grid( (void *) buf, zdf_float32, &info, &iter, "CHARGE" );
}
