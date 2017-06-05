#include <stdlib.h>
#include <stdio.h>

#include "charge.h"
#include "zdf.h"
#include "fft.h"
#include "filter.h"



void charge_new( t_charge *charge, const unsigned nx[], t_fld box[], float dt)
{

	if (( nx[0] % 2 ) || ( nx[1] % 2 )) {
		fprintf(stderr,"(*error*) Only even grid sizes are supported.");
		exit(-1);
	}

	// Number of guard cells for linear interpolation
	const unsigned int gc[2][2] = {{1,2},
		            {1,2}};

	// The FFT result will be transposed
	const unsigned int fnx[2] = { nx[1], nx[0]/2+1 };

	// Initialize grids
	scalar_grid2d_init( &charge->rho, nx, gc );
	cscalar_grid2d_init( &charge->frho, fnx, NULL );

	// Initializ FFT transform
	fftr2d_init_cfg( &charge -> fft_forward, nx[0], nx[1],
	    charge -> rho.nrow, FFT_FORWARD );

	// Set cell sizes and box limits
	unsigned i;
	for(i = 0; i<2; i++){
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

void charge_init_neutral_bkg( t_charge *charge )
{
	scalar_grid2d_init( &charge->neutral, (unsigned int *) charge ->rho.nx,
					    charge ->rho.gc );
	scalar_grid2d_zero( &charge -> neutral );
}

void charge_delete( t_charge *charge )
{
	scalar_grid2d_cleanup( &charge -> rho );
	cscalar_grid2d_cleanup( &charge -> frho );

	if ( charge -> neutral.nx[0] > 0 ) {
		scalar_grid2d_cleanup( &charge -> neutral );
	}

	fftr2d_cleanup_cfg( &charge -> fft_forward );

}

void charge_zero( t_charge *charge )
{
	scalar_grid2d_zero( &charge -> rho );
}

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


void charge_report( const t_charge *charge )
{
	char vfname[] = "charge density";

	float *buf, *p, *f;
	unsigned int i, j;

	// Pack the information
	buf = malloc( charge->rho.nx[0]*charge->rho.nx[1]*sizeof(float) );
    p = buf;
	f = charge -> rho.s;

	for( j = 0; j < charge->rho.nx[1]; j++) {
		for ( i = 0; i < charge->rho.nx[0]; i++ ) {
			p[i] = f[i];
		}
		p += charge->rho.nx[0];
		f += charge->rho.nrow;
	}


    t_zdf_grid_axis axis[2];
    axis[0] = (t_zdf_grid_axis) {
    	.min = 0.0,
    	.max = charge->box[0],
    	.label = "x_1",
    	.units = "c/\\omega_p"
    };
    axis[1] = (t_zdf_grid_axis) {
    	.min = 0.0,
    	.max = charge->box[1],
    	.label = "x_1",
    	.units = "c/\\omega_p"
    };

    t_zdf_grid_info info = {
    	.ndims = 2,
    	.label = vfname,
    	.units = "e \\omega_p^2 / c",
    	.axis = axis
    };

    info.nx[0] = charge->rho.nx[0];
    info.nx[1] = charge->rho.nx[1];

    t_zdf_iteration iter = {
    	.n = charge->iter,
    	.t = charge -> iter * charge -> dt,
    	.time_units = "1/\\omega_p"
    };

	zdf_save_grid( buf, &info, &iter, "CHARGE" );

	// free local data
	free( buf );


}
