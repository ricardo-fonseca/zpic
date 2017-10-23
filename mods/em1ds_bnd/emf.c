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


void emf_new( t_emf *emf, int nx, t_fld box, const float dt, t_fftr_cfg *fft )
{   	

	if ( nx % 2 ) {
		fprintf(stderr,"(*error*) Only even grid sizes are supported.");
		exit(-1);
	}

	// Number of guard cells for linear interpolation
	unsigned int gc[2] = {0,1}; 

	// Store pointer to required FFT configurations
	emf -> fft = fft;

	// Allocate arrays
	vfld_grid_init( &emf->E, nx, gc );
	vfld_grid_init( &emf->B, nx, gc );

	cvfld_grid_init( &emf->fEl, nx + 1, NULL );
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
	
}

void emf_delete( t_emf *emf )
{
	vfld_grid_cleanup( &emf->E );
	vfld_grid_cleanup( &emf->B );
	
	cvfld_grid_cleanup( &emf->fEl );
	cvfld_grid_cleanup( &emf->fEt );
	cvfld_grid_cleanup( &emf->fB );
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
	unsigned int i;
	
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

	fftr_r2c( &emf->fft[FORWARD], Ey, buffer );
	for( i = 1; i < emf->fEt.nx/2; i++)
		emf->fEt.y[i] += buffer[i];

	fftr_r2c( &emf->fft[FORWARD], Ez, buffer );
	for( i = 1; i < emf->fEt.nx/2; i++)
		emf->fEt.z[i] += buffer[i];

	fftr_r2c( &emf->fft[FORWARD], By, buffer );
	for( i = 1; i < emf->fEt.nx/2; i++)
		emf->fB.y[i] += buffer[i];

	fftr_r2c( &emf->fft[FORWARD], Bz, buffer );
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
	float * buf;
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
 
 Field solver
 
 *********************************************************************************************/

void advance_fB( t_emf *emf, const float dt )
{
	unsigned int i;
	
    // float complex * const restrict fBx = emf -> fB.x;
    float complex * const restrict fBy = emf -> fB.y;
    float complex * const restrict fBz = emf -> fB.z;

    //float complex * const restrict fEtx = emf -> fEt.x;
    float complex * const restrict fEty = emf -> fEt.y;
    float complex * const restrict fEtz = emf -> fEt.z;
	
    const float dk = fft_dk( emf->E.nx, emf->dx );

	// Canonical implementation
	for (i = 0; i < emf -> fB.nx; i++) {

		float kx = i * dk;

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

    // float complex * const restrict fJtx = current -> fJt.x;
    float complex * const restrict fJty = current -> fJt.y;
    float complex * const restrict fJtz = current -> fJt.z;

    const float dk = fft_dk( emf->E.nx, emf->dx );
	unsigned int i;

	for (i=0; i < emf->fEt.nx; i++) {

		float kx = i * dk;

		// fJtx is always 0 in 1D so this is unnecessary
		// fEtx[i] += dt * (                  -  fJtx[i] );  
		fEty[i] += dt * ( -I * fBz[i] * kx -  fJty[i] );  
		fEtz[i] += dt * ( +I * fBy[i] * kx -  fJtz[i] );  
 
	}
}

void update_fEl( t_emf *emf, const t_charge *charge )
{
    float complex * const restrict frho = charge -> frho.s;

    float complex * const restrict fElx = emf -> fEl.x;
    float complex * const restrict fEly = emf -> fEl.y;
    float complex * const restrict fElz = emf -> fEl.z;
    
    // For the boundary conditions we duplicate the number of cells in the charge
    const float dk = fft_dk( 2 * emf->E.nx, emf->dx );
	unsigned int i;

	fElx[0] = 0;
	fEly[0] = 0;
	fElz[0] = 0;

    for( i = 1; i < emf -> fEl.nx; i++) {
    	float kx = i * dk;
    	fElx[i] = - I * frho[i] / kx;
    	fEly[i] = 0;
    	fElz[i] = 0;
    }

}


// This code operates with periodic boundaries
void emf_update( t_emf *emf )
{
	int i;

	// Update E field

   	fftr_c2r( &emf -> fft[BACKWARD], emf -> fEt.x, emf ->E.x );
   	fftr_c2r( &emf -> fft[BACKWARD], emf -> fEt.y, emf ->E.y );
   	fftr_c2r( &emf -> fft[BACKWARD], emf -> fEt.z, emf ->E.z );

   	// Add longitudinal component
   	float * Ex_bnd = malloc( 2 * emf -> E.nx * sizeof (float));
   	fftr_c2r( &emf -> fft[BACKWARD_2], emf -> fEl.x, Ex_bnd );
   	for( i = 0; i < emf -> E.nx ; i++ ) {
   		emf ->E.x[i] += Ex_bnd[i];
   	}
   	free( Ex_bnd );


   	// Update B field
   	fftr_c2r( &emf -> fft[BACKWARD], emf -> fB.x, emf ->B.x );
   	fftr_c2r( &emf -> fft[BACKWARD], emf -> fB.y, emf ->B.y );
   	fftr_c2r( &emf -> fft[BACKWARD], emf -> fB.z, emf ->B.z );


   	// Update guard cells

    float* const restrict Ex = emf -> E.x;
    float* const restrict Ey = emf -> E.y;
    float* const restrict Ez = emf -> E.z;

    float* const restrict Bx = emf -> B.x;
    float* const restrict By = emf -> B.y;
    float* const restrict Bz = emf -> B.z;
	
    const unsigned nx = emf -> E.nx;

	// lower
	for (i = - emf->E.gc[0]; i<0; i++) {
		Ex[ i ] = Ex[ nx + i ];
		Ey[ i ] = Ey[ nx + i ];
		Ez[ i ] = Ez[ nx + i ];

		Bx[ i ] = Bx[ nx + i ];
		By[ i ] = By[ nx + i ];
		Bz[ i ] = Bz[ nx + i ];
	}

	// upper
	for (i=0; i<emf->E.gc[1]; i++) {
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
	advance_fB( emf, dt/2.0f );

	advance_fEt( emf, current, dt );

	advance_fB( emf, dt/2.0f );

	// Calculate fEl
	update_fEl( emf, charge );

	// Update (real) E, B (also updates guard cells)
	emf_update( emf );

	// Advance internal iteration number
    emf -> iter += 1;
	
    // Update timing information
	_emf_time += timer_interval_seconds(t0, timer_ticks());
}

void emf_get_energy( const t_emf *emf, double energy[] )
{
	unsigned i;

	for( i = 0; i<6; i++) energy[i] = 0;

	for( i = 0; i < emf -> E.nx; i ++ ) {
		energy[0] += emf -> E.x[i] * emf -> E.x[i];
		energy[1] += emf -> E.y[i] * emf -> E.y[i];
		energy[2] += emf -> E.z[i] * emf -> E.z[i];
		energy[3] += emf -> B.x[i] * emf -> B.x[i];
		energy[4] += emf -> B.y[i] * emf -> B.y[i];
		energy[5] += emf -> B.z[i] * emf -> B.z[i];
	}

	for( i = 0; i<6; i++) energy[i] *= 0.5 * emf -> dx;

}
