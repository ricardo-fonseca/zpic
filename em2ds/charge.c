#include <stdlib.h>
#include <stdio.h>

#include "charge.h"
#include "zdf.h"
#include "fft.h"
#include "filter.h"

/**
 * @brief Initializes Electric charge density object
 * 
 * @param charge 	Electric charge density
 * @param nx 		Number of cells [x,y]
 * @param box 		Physical box size [x,y]
 * @param dt 		Simulation time step
 */
void charge_new( t_charge *charge, const int nx[], float box[], float dt)
{

	if (( nx[0] % 2 ) || ( nx[1] % 2 )) {
		fprintf(stderr,"(*error*) Only even grid sizes are supported.");
		exit(-1);
	}

	// Number of guard cells for linear interpolation
	const unsigned int gc[2][2] = {{1,2},
		                           {1,2}};

	// The FFT result will be transposed [ky, kx]
	const unsigned int fnx[2] = { nx[1], nx[0]/2+1 };

	// Initialize grids
	scalar_grid2d_init( &charge->rho, (unsigned int *)nx, gc );
	cscalar_grid2d_init( &charge->frho, fnx, NULL );

	// Initializ FFT transform
	fftr2d_init_cfg( &charge -> fft_forward, nx[0], nx[1],
	    charge -> rho.nrow, FFT_FORWARD );

	// Set cell sizes and box limits
	for(int i = 0; i<2; i++){
		charge -> box[i] = box[i];
		charge -> dx[i]  = box[i] / nx[i];
	}

	// Initialize time information
	charge -> iter = 0;
	charge -> dt = dt;

    // Zero initial charge
    // This is only relevant for diagnostics, charge is always zeroed before deposition
	scalar_grid2d_zero( &charge -> rho );

	// Default is not to have a neutralizing background
	charge -> neutral.nx[0] = 0;
}

/**
 * @brief Initializes neutralizing background structures
 * 
 * @param charge 	Electric charge density
 */
void charge_init_neutral_bkg( t_charge *charge )
{
	scalar_grid2d_init( &charge->neutral,
						(const unsigned int *) charge ->rho.nx,
					    (const unsigned int (*)[2]) charge ->rho.gc );
	scalar_grid2d_zero( &charge -> neutral );
}

/**
 * @brief Frees dynamic memory from electric charge density
 * 
 * @param charge 	Electric charge density
 */
void charge_delete( t_charge *charge )
{
	scalar_grid2d_cleanup( &charge -> rho );
	cscalar_grid2d_cleanup( &charge -> frho );

	if ( charge -> neutral.nx[0] > 0 ) {
		scalar_grid2d_cleanup( &charge -> neutral );
	}

	fftr2d_cleanup_cfg( &charge -> fft_forward );

}

/**
 * @brief Sets all electric charge density values to zero
 * 
 * @param charge 	Electric charge density
 */
void charge_zero( t_charge *charge )
{
	scalar_grid2d_zero( &charge -> rho );
}

/**
 * @brief Advances electric charge density 1 time step
 * 
 * The routine will:
 * 1. Update the guard cells
 * 2. Add neutralizing background (if configured)
 * 3. Get the fourier transform of the charge
 * 4. Apply spectral filtering
 * 
 * @param charge 	Electric charge density
 */
void charge_update( t_charge *charge )
{
	int i, j;
	const int nrow = charge->rho.nrow;

	float* restrict const rho = charge -> rho.s;
	float complex * restrict const frho = charge -> frho.s;

	int const nx0  = charge -> rho.nx[0];
	int const nx1  = charge -> rho.nx[1];

	int const gc00 = charge -> rho.gc[0][0];
	int const gc01 = charge -> rho.gc[0][1];
	int const gc10 = charge -> rho.gc[1][0];
	int const gc11 = charge -> rho.gc[1][1];


	// x
	for (j = -gc10; j < nx1 + gc11; j++) {

		// lower - add the values from upper boundary ( both gc and inside box )
		for (i = -gc00; i < gc01; i++) {
			rho[ i + j*nrow ] += rho[ nx0 + i + j*nrow ];
		}

		// upper - just copy the values from the lower boundary
		for (i = -gc00; i < gc01; i++) {
			rho[ nx0 + i + j*nrow ] = rho[ i + j*nrow ];
		}

	}

	// y
	for (i = -gc00; i < nx0 + gc01; i++) {

		// lower - add the values from upper boundary ( both gc and inside box )
		for (j=-gc10; j<gc11; j++) {
			rho[ i + j*nrow ] += rho[ i + (nx1+j)*nrow ];
		}

		// upper - just copy the values from the lower boundary
		for (j=-gc10; j<gc11; j++) {
			rho[ i + (nx1+j)*nrow ] = rho[ i + j*nrow ];
		}

	}

	// Add neutralizing background
	// This is preferable to initializing rho to this value before charge deposition
	// because it leads to less roundoff errors
	if ( charge -> neutral.nx[0] > 0 ) {
		scalar_grid2d_add( &charge -> rho, &charge -> neutral);
	}

	// Calculate frho
	fftr2d_r2c( &charge->fft_forward, rho, frho );

	const float cutoff[2] = {0.5f,0.5f};
	cscalar2d_r2c_filter( &charge -> frho, cutoff );

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
	int i,j;
	const int nrow = charge->neutral.nrow;

	float* restrict const neutral = charge -> neutral.s;

	int const nx0  = charge -> neutral.nx[0];
	int const nx1  = charge -> neutral.nx[1];

	int const gc00 = charge -> neutral.gc[0][0];
	int const gc01 = charge -> neutral.gc[0][1];
	int const gc10 = charge -> neutral.gc[1][0];
	int const gc11 = charge -> neutral.gc[1][0];


	// x
	for (j = -gc10; j < nx1 + gc11; j++) {

		// lower - add the values from upper boundary ( both gc and inside box )
		for (i = -gc00; i < gc01; i++) {
			neutral[ i + j*nrow ] += neutral[ nx0 + i + j*nrow ];
		}

		// upper - just copy the values from the lower boundary
		for (i = -gc00; i < gc01; i++) {
			neutral[ nx0 + i + j*nrow ] = neutral[ i + j*nrow ];
		}

	}

	// y
	for (i = -gc00; i < nx0 + gc01; i++) {

		// lower - add the values from upper boundary ( both gc and inside box )
		for (j=-gc10; j<gc11; j++) {
			neutral[ i + j*nrow ] += neutral[ i + (nx1+j)*nrow ];
		}

		// upper - just copy the values from the lower boundary
		for (j=-gc10; j<gc11; j++) {
			neutral[ i + (nx1+j)*nrow ] = neutral[ i + j*nrow ];
		}

	}
	// Change sign
	for (j = -gc10; j < nx1 + gc11; j++) {
		for (i = -gc00; i < nx0 + gc01; i++) {
			neutral[ i + j*nrow ] = -neutral[ i + j*nrow ];
		}
	}

}

/**
 * @brief Saves electric charge density diagnostic information to disk
 * 
 * Saves the charge density to disk in directory "CHARGE"
 * 
 * @param charge 	Electric charge density
 */
void charge_report( const t_charge *charge )
{
	char vfname[] = "charge density";
	char vflabel[] = "\\rho";

	float *buf, *p, *f;

	// Pack the information
	buf = malloc( charge->rho.nx[0]*charge->rho.nx[1]*sizeof(float) );
    p = buf;
	f = charge -> rho.s;

	for( int j = 0; j < charge->rho.nx[1]; j++) {
		for ( int i = 0; i < charge->rho.nx[0]; i++ ) {
			p[i] = f[i];
		}
		p += charge->rho.nx[0];
		f += charge->rho.nrow;
	}

    t_zdf_grid_axis axis[2];
    axis[0] = (t_zdf_grid_axis) {
    	.min = 0.0,
    	.max = charge->box[0],
		.name = "x",
    	.label = "x",
    	.units = "c/\\omega_p"
    };
    axis[1] = (t_zdf_grid_axis) {
    	.min = 0.0,
    	.max = charge->box[1],
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

    info.count[0] = charge->rho.nx[0];
    info.count[1] = charge->rho.nx[1];

    t_zdf_iteration iter = {
		.name = "ITERATION",
    	.n = charge->iter,
    	.t = charge -> iter * charge -> dt,
    	.time_units = "1/\\omega_p"
    };

	zdf_save_grid( (void *) buf, zdf_float32, &info, &iter, "CHARGE" );

	// free local data
	free( buf );
}
