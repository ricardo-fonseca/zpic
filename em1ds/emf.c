/*
 *  emf.c
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

#include "emf.h"
#include "zdf.h"
#include "timer.h"


static double _emf_time = 0.0;

double emf_time()
{
	return _emf_time;
}

/*********************************************************************************************

 Constructor / Destructor

 *********************************************************************************************/


void emf_new( t_emf *emf, int nx, t_fld box, const float dt, t_fftr_cfg *fft_forward,
	t_fftr_cfg *fft_backward, t_filter *filter )
{

	if ( nx % 2 ) {
		fprintf(stderr,"(*error*) Only even grid sizes are supported.");
		exit(-1);
	}

	// Number of guard cells for linear interpolation
	int gc[2] = {0,1};

	// Store pointer to required FFT configuration
	emf -> fft_forward = fft_forward;
	emf -> fft_backward = fft_backward;

	// Store pointer to spectral filter data
	emf -> filter = filter;

	// Allocate arrays
	vfld_grid_init( &emf->E, nx, gc );
	vfld_grid_init( &emf->B, nx, gc );

	cvfld_grid_init( &emf->fEl, nx/2+1, NULL );
	cvfld_grid_init( &emf->fEt, nx/2+1, NULL );
	cvfld_grid_init( &emf->fB,  nx/2+1, NULL );

	// zero fields
	vfld_grid_zero( &emf->E );
	vfld_grid_zero( &emf->B );

	cvfld_grid_zero( &emf->fEl );
	cvfld_grid_zero( &emf->fEt );
	cvfld_grid_zero( &emf->fB  );

	// Set cell sizes and box limits
	emf -> box = box;
	emf -> dx = box / nx;

	// Set time step
	emf -> dt = dt;

	// Reset iteration number
	emf -> iter = 0;

	// Disable external fields by default
	emf -> ext_fld.type = EMF_EXT_FLD_NONE;
	emf -> E_part = &emf->E;
	emf -> B_part = &emf->B;

	// Set default solver type to Pseudo-spectral Analytic Time Domain
	emf -> solver_type = EMF_SOLVER_PSATD;


}

void emf_delete( t_emf *emf )
{
	vfld_grid_cleanup( &emf->E );
	vfld_grid_cleanup( &emf->B );

	cvfld_grid_cleanup( &emf->fEl );
	cvfld_grid_cleanup( &emf->fEt );
	cvfld_grid_cleanup( &emf->fB );

	if ( emf -> ext_fld.type > EMF_EXT_FLD_NONE ) {
		vfld_grid_cleanup( &emf->ext_fld.E_part_buf );
		vfld_grid_cleanup( &emf->ext_fld.B_part_buf );
		emf->E_part = NULL;
		emf->B_part = NULL;
	}

	emf -> filter = NULL;

}


/*********************************************************************************************

 Laser Pulses

*********************************************************************************************/


t_fld lon_env( const t_emf_laser* const laser, const t_fld z )
{
	if ( z > -2*laser->fwhm && z < 2*laser->fwhm ) {
		t_fld e = cos( M_PI_4 * z / laser->fwhm );
		return e*e;
	} else {
		return 0.0;
	}
}

void emf_add_laser( t_emf* const emf, const t_emf_laser* const laser )
{
	int i;

	t_fld z_center, z;
	t_fld amp, lenv, k;
	t_fld dx;
	t_fld cos_pol, sin_pol;

	float* Ey = (float *) malloc( emf->E.nx * sizeof(float));
	float* Ez = (float *) malloc( emf->E.nx * sizeof(float));
	float* By = (float *) malloc( emf->E.nx * sizeof(float));
	float* Bz = (float *) malloc( emf->E.nx * sizeof(float));

  float complex *buffer = (float complex *) malloc( emf->fEt.nx * sizeof(float complex));

	dx = emf -> dx;

	z_center = laser->start - laser->fwhm/2;
	amp = laser->omega0 * laser->a0;

	cos_pol = cos( laser -> polarization );
	sin_pol = sin( laser -> polarization );

	k = laser -> omega0;

	for (i = 0; i < emf->E.nx; i++) {

		z = i * dx - z_center;

		lenv   = amp*lon_env( laser, z );

		// Ex[i] = 0.0
		Ey[i] = +lenv * cos( k * z ) * cos_pol;
		Ez[i] = +lenv * cos( k * z ) * sin_pol;

		// Bx[i] = 0.0
		By[i] = -lenv * cos( k * z ) * sin_pol;
		Bz[i] = +lenv * cos( k * z ) * cos_pol;

		// Add laser field to real fields for diagnostic purposes
		emf->E.y[i] += Ey[i];
		emf->E.z[i] += Ez[i];
		emf->B.y[i] += By[i];
		emf->B.z[i] += Bz[i];
	}

	// Filter laser field by only adding components up to
	// half the Nyquist limit

	fftr_r2c( emf->fft_forward, Ey, buffer );
	for( i = 1; i < emf->fEt.nx/2; i++)
		emf->fEt.y[i] += buffer[i];

	fftr_r2c( emf->fft_forward, Ez, buffer );
	for( i = 1; i < emf->fEt.nx/2; i++)
		emf->fEt.z[i] += buffer[i];

	fftr_r2c( emf->fft_forward, By, buffer );
	for( i = 1; i < emf->fEt.nx/2; i++)
		emf->fB.y[i] += buffer[i];

	fftr_r2c( emf->fft_forward, Bz, buffer );
	for( i = 1; i < emf->fEt.nx/2; i++)
		emf->fB.z[i] += buffer[i];

	free( buffer );
	free( Bz );
	free( By );
	free( Ez );
	free( Ey );

}

/*********************************************************************************************

 Diagnostics

 *********************************************************************************************/


void emf_report( const t_emf *emf, const char field, const char fc )
{
	char vfname[3];

	// Choose field to save
	const t_vfld_grid * f;
	switch (field) {
		case EFLD:
			f = &emf->E;
			vfname[0] = 'E';
			break;
		case BFLD:
			f = &emf->B;
			vfname[0] = 'B';
			break;
		default:
			printf("Invalid field type selected, returning\n");
			return;
	}

	// Pack the information
	float * buf = NULL;
	switch (fc) {
		case 0:
			buf = f->x;
			vfname[1] = '1';
			break;
		case 1:
			buf = f->y;
			vfname[1] = '2';
			break;
		case 2:
			buf = f->z;
			vfname[1] = '3';
			break;
	}
	vfname[2] = 0;

    t_zdf_grid_axis axis[1];
    axis[0] = (t_zdf_grid_axis) {
    	.min = 0.0,
    	.max = emf->box,
    	.label = "x_1",
    	.units = "c/\\omega_p"
    };

    t_zdf_grid_info info = {
    	.ndims = 1,
    	.label = vfname,
    	.units = "m_e c \\omega_p e^{-1}",
    	.axis = axis
    };

    info.nx[0] = emf->E.nx;

    t_zdf_iteration iter = {
    	.n = emf->iter,
    	.t = emf -> iter * emf -> dt,
    	.time_units = "1/\\omega_p"
    };

	zdf_save_grid( buf, &info, &iter, "EMF" );

}


/*********************************************************************************************

 Spectral Field solver

 *********************************************************************************************/

void advance_fB( t_emf *emf, const float dt )
{
	int i;

    // float complex * const restrict fBx = emf -> fB.x;
    float complex * const restrict fBy = emf -> fB.y;
    float complex * const restrict fBz = emf -> fB.z;

    //float complex * const restrict fEtx = emf -> fEt.x;
    float complex * const restrict fEty = emf -> fEt.y;
    float complex * const restrict fEtz = emf -> fEt.z;

    const float dk = fft_dk( emf->E.nx, emf->dx );

	// Canonical implementation
	for (i = 0; i < emf -> fB.nx; i++) {

		const float kx = i * dk;

		// fBx[i] += 0;
		fBy[i] +=  + dt * I * fEtz[i] * kx;
		fBz[i] +=  - dt * I * fEty[i] * kx;
	}
}

void advance_fEt( t_emf *emf, const t_current *current, const float dt )
{

    // float complex * const restrict fEtx = emf -> fEt.x;
    float complex * const restrict fEty = emf -> fEt.y;
    float complex * const restrict fEtz = emf -> fEt.z;

    // float complex * const restrict fBx = emf -> fB.x;
    float complex * const restrict fBy = emf -> fB.y;
    float complex * const restrict fBz = emf -> fB.z;

    // float complex * const restrict fJx = current -> fJ.x;
    float complex * const restrict fJy = current -> fJ.y;
    float complex * const restrict fJz = current -> fJ.z;

    const float dk = fft_dk( emf->E.nx, emf->dx );

	// kx = 0
	{
		// const float complex fJtx = 0;
		const float complex fJty = fJy[0];
		const float complex fJtz = fJz[0];

		// fEtx[0] -= dt * fJtx;
		fEty[0] -= dt * fJty;
		fEtz[0] -= dt * fJtz;
	}

	// kx > 1
	for (int i=1; i < emf->fEt.nx; i++) {

		const float kx = i * dk;

		// Get transverse current
		//const float complex fJtx = 0;
		const float complex fJty = fJy[i];
		const float complex fJtz = fJz[i];

		// Advance field
		//fEtx[i] += dt * (                  -  fJtx );
		fEty[i] += dt * ( -I * fBz[i] * kx -  fJty );
		fEtz[i] += dt * ( +I * fBy[i] * kx -  fJtz );

	}
}

void update_fEl( t_emf *emf, const t_charge *charge )
{
    float complex * const restrict frho = charge -> frho.s;

    float complex * const restrict fElx = emf -> fEl.x;
    float complex * const restrict fEly = emf -> fEl.y;
    float complex * const restrict fElz = emf -> fEl.z;

    const float dk = fft_dk( emf->E.nx, emf->dx );
	int i;

	// kx = 0
	fElx[0] = 0;
	fEly[0] = 0;
	fElz[0] = 0;

	// kx > 0
    for( i = 1; i < emf -> fEl.nx; i++) {
    	const float kx = i * dk;
    	fElx[i] = - I * frho[i] / kx;
    	fEly[i] = 0;
    	fElz[i] = 0;
    }

}

/*********************************************************************************************

 Pseudo-spectral Analytical time domain field solver

 *********************************************************************************************/

void advance_psatd( t_emf *emf, const t_current *current, const float dt )
{

    // float complex * const restrict fEtx = emf -> fEt.x;
    float complex * const restrict fEty = emf -> fEt.y;
    float complex * const restrict fEtz = emf -> fEt.z;

    // float complex * const restrict fBx = emf -> fB.x;
    float complex * const restrict fBy = emf -> fB.y;
    float complex * const restrict fBz = emf -> fB.z;

    // float complex * const restrict fJtx = current -> fJt.x;
    float complex * const restrict fJy = current -> fJ.y;
    float complex * const restrict fJz = current -> fJ.z;

    const float dk = fft_dk( emf->E.nx, emf->dx );

    // kx = 0
    {
		//const float complex fJtx = 0;
		const float complex fJty = fJy[0];
		const float complex fJtz = fJz[0];

		//fEtx[0] -= dt * fJtx;
		fEty[0] -= dt * fJty;
		fEtz[0] -= dt * fJtz;
    }

	// kx > 0
	for (int i=1; i < emf->fEt.nx; i++) {

		const float kx = i * dk;

		// Get transverse current
		//const float complex fJtx = 0;
		const float complex fJty = fJy[i];
		const float complex fJtz = fJz[i];

		// Advance field

		const float C = cosf( kx * dt );
		const float S = sinf( kx * dt );

		float complex Ey = fEty[i];
		float complex Ez = fEtz[i];

		float complex By = fBy[i];
		float complex Bz = fBz[i];

		// fJtx is always 0 in 1D so this is unnecessary
		// Ex = C * Ex;
		Ey = C * Ey - I * S * fBz[i] - S * fJty[i] / kx;
		Ez = C * Ez + I * S * fBy[i] - S * fJtz[i] / kx;

		// Bx = C * Bx;
		By = C * By + I * S * fEtz[i] - I * (1.0f - C) * fJtz[i] / kx;
		Bz = C * Bz - I * S * fEty[i] + I * (1.0f - C) * fJty[i] / kx;

		// fEtx[i] = Ex;
		fEty[i] = Ey;
        fEtz[i] = Ez;

        // fBx[i]  = Bx;
        fBy[i]  = By;
        fBz[i]  = Bz;

	}
}


// This code operates with periodic boundaries
void emf_update( t_emf *emf )
{
	int i;

	// Update E field
	// To save memory we store the complete E field (longitudinal + transverse)
	// in the longitudinal E field array

	if ( emf -> filter -> type > FILTER_NONE ) {
	    float * const restrict Sk = emf -> filter -> Sk;

		for( i = 0; i < emf -> fEl.nx; i++) {
			emf -> fEl.x[i] = Sk[i] * (emf ->fEl.x[i] + emf -> fEt.x[i]);
			emf -> fEl.y[i] = Sk[i] * (emf ->fEl.y[i] + emf -> fEt.y[i]);
			emf -> fEl.z[i] = Sk[i] * (emf ->fEl.z[i] + emf -> fEt.z[i]);
		}
	} else {
		for( i = 0; i < emf -> fEl.nx; i++) {
			emf -> fEl.x[i] = emf ->fEl.x[i] + emf -> fEt.x[i];
			emf -> fEl.y[i] = emf ->fEl.y[i] + emf -> fEt.y[i];
			emf -> fEl.z[i] = emf ->fEl.z[i] + emf -> fEt.z[i];
		}
	}

   	fftr_c2r( emf -> fft_backward, emf -> fEl.x, emf ->E.x );
   	fftr_c2r( emf -> fft_backward, emf -> fEl.y, emf ->E.y );
   	fftr_c2r( emf -> fft_backward, emf -> fEl.z, emf ->E.z );

   	// Update B field
   	fftr_c2r( emf -> fft_backward, emf -> fB.x, emf ->B.x );
   	fftr_c2r( emf -> fft_backward, emf -> fB.y, emf ->B.y );
   	fftr_c2r( emf -> fft_backward, emf -> fB.z, emf ->B.z );


   	// Update guard cells
    float* const restrict Ex = emf -> E.x;
    float* const restrict Ey = emf -> E.y;
    float* const restrict Ez = emf -> E.z;

    float* const restrict Bx = emf -> B.x;
    float* const restrict By = emf -> B.y;
    float* const restrict Bz = emf -> B.z;

    const int nx = emf -> E.nx;
    const int gc0 = emf->E.gc[0];
    const int gc1 = emf->E.gc[1];

	// lower
	for (i = -gc0; i<0; i++) {
		Ex[ i ] = Ex[ nx + i ];
		Ey[ i ] = Ey[ nx + i ];
		Ez[ i ] = Ez[ nx + i ];

		Bx[ i ] = Bx[ nx + i ];
		By[ i ] = By[ nx + i ];
		Bz[ i ] = Bz[ nx + i ];
	}

	// upper
	for (i=0; i<gc1; i++) {
		Ex[ nx + i ] = Ex[ i ];
		Ey[ nx + i ] = Ey[ i ];
		Ez[ nx + i ] = Ez[ i ];

		Bx[ nx + i ] = Bx[ i ];
		By[ nx + i ] = By[ i ];
		Bz[ nx + i ] = Bz[ i ];
	}

}

void emf_advance( t_emf *emf, const t_charge *charge, const t_current *current )
{
	uint64_t t0 = timer_ticks();
	const float dt = emf->dt;

	// Advance fB, fEt

	if ( emf -> solver_type == EMF_SOLVER_PSATD ) {
		// Pseudo-spectral Analytical Time-Domain
		advance_psatd( emf, current, dt );
	} else {
		// Pseudo-spectral Time Domain
		advance_fB( emf, dt/2.0f );
		advance_fEt( emf, current, dt );
		advance_fB( emf, dt/2.0f );
	}

	// Calculate fEl
	update_fEl( emf, charge );

	// Update (real) E, B (also updates guard cells)
	emf_update( emf );

	// Update contribuition of external fields if necessary
	if ( emf -> ext_fld.type > EMF_EXT_FLD_NONE ) emf_update_part_fld( emf );

	// Advance internal iteration number
    emf -> iter += 1;

    // Update timing information
	_emf_time += timer_interval_seconds(t0, timer_ticks());
}

void emf_get_energy( const t_emf *emf, double energy[] )
{
	int i;
    float* const restrict Ex = emf -> E.x;
    float* const restrict Ey = emf -> E.y;
    float* const restrict Ez = emf -> E.z;

    float* const restrict Bx = emf -> B.x;
    float* const restrict By = emf -> B.y;
    float* const restrict Bz = emf -> B.z;

	for( i = 0; i<6; i++) energy[i] = 0;

	for( i = 0; i < emf -> E.nx; i ++ ) {
		energy[0] += Ex[i] * Ex[i];
		energy[1] += Ey[i] * Ey[i];
		energy[2] += Ez[i] * Ez[i];
		energy[3] += Bx[i] * Bx[i];
		energy[4] += By[i] * By[i];
		energy[5] += Bz[i] * Bz[i];
	}

	for( i = 0; i<6; i++) energy[i] *= 0.5 * emf -> dx;

}

/*********************************************************************************************

External Fields

 *********************************************************************************************/

int emf_set_ext_fld( t_emf* const emf, t_emf_ext_fld* ext_fld ) {

	emf -> ext_fld.type = ext_fld -> type;

	if ( emf -> ext_fld.type == EMF_EXT_FLD_NONE ) {
		// Particle fields just point to the self-consistent fields
		emf -> E_part = &(emf -> E);
		emf -> B_part = &(emf -> B);
	} else {
	    switch( emf -> ext_fld.type ) {
	        case( EMF_EXT_FLD_UNIFORM ):
	        	emf -> ext_fld.E0 = ext_fld->E0;
	        	emf -> ext_fld.B0 = ext_fld->B0;
	        	break;

	    	default:
	    		fprintf(stderr, "Invalid external field type, aborting.\n" );
				return -1;
	    }

		// Allocate space for additional field grids
		vfld_grid_init( &emf -> ext_fld.E_part_buf, emf -> E.nx, emf -> E.gc );
		vfld_grid_init( &emf -> ext_fld.B_part_buf, emf -> B.nx, emf -> B.gc );

		emf -> E_part = &emf->ext_fld.E_part_buf;
	    emf -> B_part = &emf->ext_fld.B_part_buf;

	    // Initialize values on E/B_part grids
	    emf_update_part_fld( emf );

	}

	return 0;

}

/**
 * Updates field values seen by particles with externally imposed fields
 * @param emf EMF object holding field data
 */
void emf_update_part_fld( t_emf* const emf ) {

    float* const restrict E_part_x = emf -> ext_fld.E_part_buf.x;
    float* const restrict E_part_y = emf -> ext_fld.E_part_buf.y;
    float* const restrict E_part_z = emf -> ext_fld.E_part_buf.z;

    float* const restrict B_part_x = emf -> ext_fld.B_part_buf.x;
    float* const restrict B_part_y = emf -> ext_fld.B_part_buf.y;
    float* const restrict B_part_z = emf -> ext_fld.B_part_buf.z;

    // Currently only EMF_EXT_FLD_UNIFORM is supported

	// Add external field values to self consistent fields
	for( int i = -emf->E.gc[0]; i < emf->E.nx + emf->E.gc[1]; i++ ){
		E_part_x[i] = emf -> E.x[i] + emf -> ext_fld.E0.x;
		E_part_y[i] = emf -> E.y[i] + emf -> ext_fld.E0.y;
		E_part_z[i] = emf -> E.z[i] + emf -> ext_fld.E0.z;

		B_part_x[i] = emf -> B.x[i] + emf -> ext_fld.B0.x;
		B_part_y[i] = emf -> B.y[i] + emf -> ext_fld.B0.y;
		B_part_z[i] = emf -> B.z[i] + emf -> ext_fld.B0.z;
	}
}

