/**
 * @file emf.c
 * @author Ricardo Fonseca
 * @brief EM Fields
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

#include "emf.h"
#include "zdf.h"
#include "timer.h"

void emf_update_part_fld( t_emf* const emf );

/// Time spent advancing the EM fields
static double _emf_time = 0.0;

/**
 * @brief Time spent advancing the EM fields
 * 
 * @return      Time spent in seconds
 */
double emf_time()
{
    return _emf_time;
}

/*********************************************************************************************

 Constructor / Destructor

 *********************************************************************************************/

/**
 * @brief Initializes the EM field object
 * 
 * @param emf			EM fields object
 * @param nx			Number of cells
 * @param box			Physical box size
 * @param dt			Simulation time step
 * @param fft_forward	FFT configuration for forward transforms (shared
 * 						with other objects)
 * @param fft_backward	FFT configuration for backward transforms (shared
 * 						with other objects)
 * @param filter		Spectral filtering parameters
 */
void emf_new( t_emf *emf, int nx, float box, const float dt, t_fftr_cfg *fft_forward,
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
    float3_grid_init( &emf->E, nx, gc );
    float3_grid_init( &emf->B, nx, gc );

    cfloat3_grid_init( &emf->fEl, nx/2+1, NULL );
    cfloat3_grid_init( &emf->fEt, nx/2+1, NULL );
    cfloat3_grid_init( &emf->fB,  nx/2+1, NULL );

    // zero fields
    float3_grid_zero( &emf->E );
    float3_grid_zero( &emf->B );

    cfloat3_grid_zero( &emf->fEl );
    cfloat3_grid_zero( &emf->fEt );
    cfloat3_grid_zero( &emf->fB  );

    // Set cell sizes and box limits
    emf -> box = box;
    emf -> dx = box / nx;

    // Set time step
    emf -> dt = dt;

    // Reset iteration number
    emf -> iter = 0;

    // Disable external fields by default
    emf -> ext_fld.E_type = EMF_FLD_TYPE_NONE;
    emf -> ext_fld.B_type = EMF_FLD_TYPE_NONE;
    emf -> E_part = &emf->E;
    emf -> B_part = &emf->B;

    // Set default solver type to Pseudo-spectral Analytic Time Domain
    emf -> solver_type = EMF_SOLVER_PSATD;


}

/**
 * @brief Frees dynamic memory from EM fields.
 * 
 * If external fields are in use, the dynamic memory associated with these
 * will also be freed.
 * 
 * @param emf 	EM fields
 */
void emf_delete( t_emf *emf )
{
    float3_grid_cleanup( &emf->E );
    float3_grid_cleanup( &emf->B );

    cfloat3_grid_cleanup( &emf->fEl );
    cfloat3_grid_cleanup( &emf->fEt );
    cfloat3_grid_cleanup( &emf->fB );

    if ( emf -> ext_fld.E_type > EMF_FLD_TYPE_NONE ) {
        float3_grid_cleanup( &emf->ext_fld.E_part_buf );
    }

    if ( emf -> ext_fld.B_type > EMF_FLD_TYPE_NONE ) {
        float3_grid_cleanup( &emf->ext_fld.B_part_buf );
    }

    emf -> filter = NULL;

}


/*********************************************************************************************

 Laser Pulses

*********************************************************************************************/


/**
 * @brief Determines longitudinal envelope value of laser pulse
 * 
 * @param laser 	Laser pulse parameters
 * @param z 		Longitudinal position
 * @return lon_env	Envelope value
 */
float lon_env( const t_emf_laser* const laser, const float z )
{
	if ( z > laser -> start ) {
		// Ahead of laser
		return 0.0;
	} else if ( z > laser -> start - laser -> rise ) {
		// Laser rise
		float csi = z - laser -> start;
		float e = sin( M_PI_2 * csi / laser->rise );
		return e*e;
	} else if ( z > laser -> start - (laser -> rise + laser -> flat) ) {
		// Flat-top
		return 1.0;
	} else if ( z > laser -> start - (laser -> rise + laser -> flat + laser -> fall) ) {
		// Laser fall
		float csi = z - (laser -> start - laser -> rise - laser -> flat - laser -> fall);
		float e = sin( M_PI_2 * csi / laser->fall );
		return e*e;
	}

	// Before laser
	return 0.0;
}

/**
 * @brief Add laser pulse to simulation.
 * 
 * Laser pulses are superimposed on top of existing E and B fields. 
 * Multiple lasers can be added.
 * 
 * @param emf 		EM fields
 * @param laser 	Laser pulse parameters
 */
void emf_add_laser( t_emf* const emf, t_emf_laser* const laser )
{

	// Validate laser parameters
	if ( laser -> fwhm != 0 ) {
		if ( laser -> fwhm <= 0 ) {
			fprintf(stderr, "Invalid laser FWHM, must be > 0, aborting.\n" );
			exit(-1);
		}

		// The fwhm parameter overrides the rise/flat/fall parameters
		laser -> rise = laser -> fwhm;
		laser -> fall = laser -> fwhm;
		laser -> flat = 0.;
	}

	if ( laser -> rise <= 0 ) {
		fprintf(stderr, "Invalid laser RISE, must be > 0, aborting.\n" );
		exit(-1);
	}

	if ( laser -> flat < 0 ) {
		fprintf(stderr, "Invalid laser FLAT, must be >= 0, aborting.\n" );
		exit(-1);
	}

	if ( laser -> fall <= 0 ) {
		fprintf(stderr, "Invalid laser FALL, must be > 0, aborting.\n" );
		exit(-1);
	}

	// Launch laser
    float z_center, z;
    float amp, lenv, k;
    float dx;
    float cos_pol, sin_pol;

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

    for (int i = 0; i < emf->E.nx; i++) {

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
    for( int i = 1; i < emf->fEt.nx/2; i++)
        emf->fEt.y[i] += buffer[i];

    fftr_r2c( emf->fft_forward, Ez, buffer );
    for( int i = 1; i < emf->fEt.nx/2; i++)
        emf->fEt.z[i] += buffer[i];

    fftr_r2c( emf->fft_forward, By, buffer );
    for( int i = 1; i < emf->fEt.nx/2; i++)
        emf->fB.y[i] += buffer[i];

    fftr_r2c( emf->fft_forward, Bz, buffer );
    for( int i = 1; i < emf->fEt.nx/2; i++)
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

/**
 * @brief Saves EM fields diagnostic information to disk
 * 
 * Saves the selected type / density component to disk in directory
 * "EMF". Guard cell values are discarded.
 *
 * @param emf 		EM Fields
 * @param field 	Which field to save (E, B, Epart, Bpart)
 * @param fc 		Field component to save, must be one of {0,1,2}
 */
void emf_report( const t_emf *emf, const char field, const int fc )
{
	char vfname[16];	// Dataset name
	char vflabel[16];	// Dataset label (for plots)

	char comp[] = {'x','y','z'};

	if ( fc < 0 || fc > 2 ) {
		fprintf(stderr, "(*error*) Invalid field component (fc) selected, returning\n");
		return;
	}

    // Choose field to save
    t_float3_grid * f;
    switch (field) {
        case EFLD:
            f = (t_float3_grid *) &emf->E;
            snprintf(vfname,16,"E%1d",fc);
            snprintf(vflabel,16,"E_%c",comp[fc]);
            break;
        case BFLD:
            f = (t_float3_grid *) &emf->B;
            snprintf(vfname,16,"E%1d",fc);
            snprintf(vflabel,16,"E_%c",comp[fc]);
            break;
		case EPART:
			f = (t_float3_grid *) &emf->E_part;
            snprintf(vfname,16,"E%1d-part",fc);
            snprintf(vflabel,16,"E_{%cp}",comp[fc]);
			break;
		case BPART:
			f = (t_float3_grid *) &emf->B_part;
            snprintf(vfname,16,"B%1d-part",fc);
            snprintf(vflabel,16,"B_{%cp}",comp[fc]);
			break;
        default:
            printf("Invalid field type selected, returning\n");
            return;
    }

    // Choose the field component
    float * buf = NULL;
    switch (fc) {
        case 0:
            buf = f->x;
            break;
        case 1:
            buf = f->y;
            break;
        case 2:
            buf = f->z;
            break;
    }

    t_zdf_grid_axis axis[1];
    axis[0] = (t_zdf_grid_axis) {
        .min = 0.0,
        .max = emf->box,
        .name = "x",
        .label = "x",
        .units = "c/\\omega_p"
    };

    t_zdf_grid_info info = {
        .ndims = 1,
		.name = vfname,
    	.label = vflabel,
        .units = "m_e c \\omega_p e^{-1}",
        .axis = axis
    };

    info.count[0] = emf->E.nx;

    t_zdf_iteration iter = {
        .name = "ITERATION",
        .n = emf->iter,
        .t = emf -> iter * emf -> dt,
        .time_units = "1/\\omega_p"
    };

    zdf_save_grid( (float *) buf, zdf_float32, &info, &iter, "EMF" );
}


/*********************************************************************************************

 Spectral Field solver

 *********************************************************************************************/

/**
 * @brief Advance the magnetic field using a PSTD algorithm
 * 
 * @param emf   EM fields
 * @param dt    Time step
 */
void advance_fB( t_emf *emf, const float dt )
{

    // float complex * const restrict fBx = emf -> fB.x;
    float complex * const restrict fBy = emf -> fB.y;
    float complex * const restrict fBz = emf -> fB.z;

    //float complex * const restrict fEtx = emf -> fEt.x;
    float complex * const restrict fEty = emf -> fEt.y;
    float complex * const restrict fEtz = emf -> fEt.z;

    const float dk = fft_dk( emf->E.nx, emf->dx );

    // Canonical implementation
    for (int i = 0; i < emf -> fB.nx; i++) {

        const float kx = i * dk;

        // fBx[i] += 0;
        fBy[i] +=  + dt * I * fEtz[i] * kx;
        fBz[i] +=  - dt * I * fEty[i] * kx;
    }
}

/**
 * @brief Advance the transverse component of the Electric field using a PSTD algorithm
 * 
 * @param emf       EM fields
 * @param current   Electric current density
 * @param dt        Time step
 */
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

/**
 * @brief Update longitudinal Electric field from charge density
 * 
 * Note: We chose to set the k=0 components to 0
 * 
 * @param emf       EM Fields
 * @param charge    Charge density
 */
void update_fEl( t_emf *emf, const t_charge *charge )
{
    float complex * const restrict frho = charge -> frho.s;

    float complex * const restrict fElx = emf -> fEl.x;
    float complex * const restrict fEly = emf -> fEl.y;
    float complex * const restrict fElz = emf -> fEl.z;

    const float dk = fft_dk( emf->E.nx, emf->dx );

    // kx = 0
    fElx[0] = 0;
    fEly[0] = 0;
    fElz[0] = 0;

    // kx > 0
    for( int i = 1; i < emf -> fEl.nx; i++) {
        const float kx = i * dk;
        fElx[i] = - I * frho[i] / kx;
        fEly[i] = 0;
        fElz[i] = 0;
    }

}

/*********************************************************************************************

 Pseudo-spectral Analytical time domain field solver

 *********************************************************************************************/

/**
 * @brief Advance transverse E and B using the PSATD algorithm
 * 
 * @param emf       EM fields
 * @param current   Electric current density
 * @param dt        Time step
 */
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
        Ey = C * Ey - I * S * fBz[i] - S * fJty / kx;
        Ez = C * Ez + I * S * fBy[i] - S * fJtz / kx;

        // Bx = C * Bx;
        By = C * By + I * S * fEtz[i] - I * (1.0f - C) * fJtz / kx;
        Bz = C * Bz - I * S * fEty[i] + I * (1.0f - C) * fJty / kx;

        // fEtx[i] = Ex;
        fEty[i] = Ey;
        fEtz[i] = Ez;

        // fBx[i]  = Bx;
        fBy[i]  = By;
        fBz[i]  = Bz;

    }
}

/**
 * @brief Updates real electric and magnetic field from Fourier transforms
 * 
 * The routine will apply spectral filtering, if supplied, and update guard
 * cell values. Note that this code operates with periodic boundaries.
 * 
 * It will also update the guard cell values for field interpolation
 * 
 * @param emf   EM fields
 */
void emf_update( t_emf *emf )
{
    // Update E field
    // To save memory we store the complete E field (longitudinal + transverse)
    // in the longitudinal E field array

    if ( emf -> filter -> type > FILTER_NONE ) {
        float * const restrict Sk = emf -> filter -> Sk;

        for( int i = 0; i < emf -> fEl.nx; i++) {
            emf -> fEl.x[i] = Sk[i] * (emf ->fEl.x[i] + emf -> fEt.x[i]);
            emf -> fEl.y[i] = Sk[i] * (emf ->fEl.y[i] + emf -> fEt.y[i]);
            emf -> fEl.z[i] = Sk[i] * (emf ->fEl.z[i] + emf -> fEt.z[i]);
        }
    } else {
        for( int i = 0; i < emf -> fEl.nx; i++) {
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
    for (int i = -gc0; i<0; i++) {
        Ex[ i ] = Ex[ nx + i ];
        Ey[ i ] = Ey[ nx + i ];
        Ez[ i ] = Ez[ nx + i ];

        Bx[ i ] = Bx[ nx + i ];
        By[ i ] = By[ nx + i ];
        Bz[ i ] = Bz[ nx + i ];
    }

    // upper
    for (int i=0; i<gc1; i++) {
        Ex[ nx + i ] = Ex[ i ];
        Ey[ nx + i ] = Ey[ i ];
        Ez[ nx + i ] = Ez[ i ];

        Bx[ nx + i ] = Bx[ i ];
        By[ nx + i ] = By[ i ];
        Bz[ nx + i ] = Bz[ i ];
    }

}

/**
 * @brief Advance EM fields 1 timestep
 * 
 * Fields are advanced in time using a spectral algorithm. The routine will also:
 * 1. Update guard cell values / apply boundary conditions
 * 2. Apply spectral filtering, if configured
 * 3. Update "particle" fields if using external fields
 * 
 * @param emf 		EM fields
 * @param charge    Electric charge density
 * @param current 	Electric current density
 */
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
    emf_update_part_fld( emf );

    // Advance internal iteration number
    emf -> iter += 1;

    // Update timing information
    _emf_time += timer_interval_seconds(t0, timer_ticks());
}

/**
 * @brief Calculate total EM field energy
 *
 * @param[in] emf EM field
 * @param[out] energy Energy values vector
 */
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
 * External Fields
 *********************************************************************************************/

/**
 * Sets the external fields to be used for the simulation
 * @param   emf     EM field object
 * @param   ext_fld External fields
 */
void emf_set_ext_fld( t_emf* const emf, t_emf_ext_fld* ext_fld ) {

    emf -> ext_fld.E_type = ext_fld -> E_type;

    if ( emf -> ext_fld.E_type == EMF_FLD_TYPE_NONE ) {
        // Particle fields just point to the self-consistent fields
        emf -> E_part = &(emf -> E);
    } else {
        switch( emf -> ext_fld.E_type ) {
            case( EMF_FLD_TYPE_UNIFORM ):
                emf -> ext_fld.E_0 = ext_fld->E_0;
                break;

            case( EMF_FLD_TYPE_CUSTOM ):
                emf -> ext_fld.E_custom = ext_fld->E_custom;
                emf -> ext_fld.E_custom_data = ext_fld->E_custom_data;
                break;

            default:
                fprintf(stderr, "Invalid external field type, aborting.\n" );
                exit(-1);
        }

        // Allocate space for additional field grids
        float3_grid_init( &emf -> ext_fld.E_part_buf, emf -> E.nx, emf -> E.gc );
        emf->E_part = &emf->ext_fld.E_part_buf;
    }

    emf -> ext_fld.B_type = ext_fld -> B_type;

    if ( emf -> ext_fld.B_type == EMF_FLD_TYPE_NONE ) {
        // Particle fields just point to the self-consistent fields
        emf -> B_part = &(emf -> B);
    } else {
        switch( emf -> ext_fld.B_type ) {
            case( EMF_FLD_TYPE_UNIFORM ):
                emf -> ext_fld.B_0 = ext_fld->B_0;
                break;

            case( EMF_FLD_TYPE_CUSTOM ):
                emf -> ext_fld.B_custom = ext_fld->B_custom;
                emf -> ext_fld.B_custom_data = ext_fld->B_custom_data;
                break;

            default:
                fprintf(stderr, "Invalid external field type, aborting.\n" );
                exit(-1);
        }

        // Allocate space for additional field grids
        float3_grid_init( &emf -> ext_fld.B_part_buf, emf -> B.nx, emf -> B.gc );
        emf -> B_part = &emf->ext_fld.B_part_buf;
    }

    // Initialize values on E/B_part grids
    emf_update_part_fld( emf );

}


/**
 * @brief Updates field values seen by particles with externally imposed fields
 * 
 * @param emf   EM field object
 */
void emf_update_part_fld( t_emf* const emf ) {

    float* const restrict E_part_x = emf -> ext_fld.E_part_buf.x;
    float* const restrict E_part_y = emf -> ext_fld.E_part_buf.y;
    float* const restrict E_part_z = emf -> ext_fld.E_part_buf.z;

    switch (emf->ext_fld.E_type)
    {
    case EMF_FLD_TYPE_UNIFORM: {
        for (int i=-emf->E.gc[0]; i<emf->E.nx+emf->E.gc[1]; i++) {
            E_part_x[i] = emf -> E.x[i] + emf -> ext_fld.E_0.x;
            E_part_y[i] = emf -> E.y[i] + emf -> ext_fld.E_0.y;
            E_part_z[i] = emf -> E.z[i] + emf -> ext_fld.E_0.z;
        }
        break; }
    case EMF_FLD_TYPE_CUSTOM: {
        for (int i=-emf->E.gc[0]; i<emf->E.nx+emf->E.gc[1]; i++) {
            float3 ext_E = (*emf->ext_fld.E_custom)(i,emf->dx,emf->ext_fld.E_custom_data);

            E_part_x[i] = emf -> E.x[i] + ext_E.x;
            E_part_y[i] = emf -> E.y[i] + ext_E.y;
            E_part_z[i] = emf -> E.z[i] + ext_E.z;
        }
        break; }
    case EMF_FLD_TYPE_NONE:
        break;
    }

    float* const restrict B_part_x = emf -> ext_fld.B_part_buf.x;
    float* const restrict B_part_y = emf -> ext_fld.B_part_buf.y;
    float* const restrict B_part_z = emf -> ext_fld.B_part_buf.z;

    switch (emf->ext_fld.B_type)
    {
    case EMF_FLD_TYPE_UNIFORM: {
        for (int i=-emf->B.gc[0]; i<emf->B.nx+emf->B.gc[1]; i++) {
            B_part_x[i] = emf -> B.x[i] + emf -> ext_fld.B_0.x;
            B_part_y[i] = emf -> B.y[i] + emf -> ext_fld.B_0.y;
            B_part_z[i] = emf -> B.z[i] + emf -> ext_fld.B_0.z;
        }
        break; }
    case EMF_FLD_TYPE_CUSTOM: {
        for (int i=-emf->B.gc[0]; i<emf->B.nx+emf->B.gc[1]; i++) {
            float3 ext_B = (*emf->ext_fld.B_custom)(i,emf->dx,emf->ext_fld.B_custom_data);

            B_part_x[i] = emf -> B.x[i] + ext_B.x;
            B_part_y[i] = emf -> B.y[i] + ext_B.y;
            B_part_z[i] = emf -> B.z[i] + ext_B.z;
        }
        break; }
    case EMF_FLD_TYPE_NONE:
        break;
    }
}

/**
 * Initialize EMF field values
 * @param emf       EM field object
 * @param init_fld  Initial field parameters
 */
void emf_init_fld( t_emf* const emf, t_emf_init_fld* init_fld )
{
    if ( emf -> iter != 0 ) {
        fprintf(stderr, "emf_init_fld should only be called at initialization, aborting...\n" );
        exit(-1);
    }

    float* const restrict Ex = emf -> E.x;
    float* const restrict Ey = emf -> E.y;
    float* const restrict Ez = emf -> E.z;

    float* const restrict Bx = emf -> B.x;
    float* const restrict By = emf -> B.y;
    float* const restrict Bz = emf -> B.z;

    switch ( init_fld -> E_type )
    {
    case EMF_FLD_TYPE_NONE:
        break;

    case EMF_FLD_TYPE_UNIFORM:
        for (int i=-emf->E.gc[0]; i<emf->E.nx+emf->E.gc[1]; i++) {
            Ex[ i ] = init_fld -> E_0.x;
            Ey[ i ] = init_fld -> E_0.y;
            Ez[ i ] = init_fld -> E_0.z;
        }

        fftr_r2c( emf->fft_forward, Ex, emf->fEt.x );
        fftr_r2c( emf->fft_forward, Ey, emf->fEt.y );
        fftr_r2c( emf->fft_forward, Ez, emf->fEt.z );

        break;

    case EMF_FLD_TYPE_CUSTOM:
        for (int i=-emf->E.gc[0]; i<emf->E.nx+emf->E.gc[1]; i++) {
            float3 init_E = (init_fld->E_custom)
                (i,emf->dx, init_fld->E_custom_data);
            Ex[ i ] = init_E.x;
            Ey[ i ] = init_E.y;
            Ez[ i ] = init_E.z;
        }

        fftr_r2c( emf->fft_forward, Ex, emf->fEt.x );
        fftr_r2c( emf->fft_forward, Ey, emf->fEt.y );
        fftr_r2c( emf->fft_forward, Ez, emf->fEt.z );

        break;
    }    

    switch ( init_fld -> B_type )
    {
    case EMF_FLD_TYPE_NONE:
        break;

    case EMF_FLD_TYPE_UNIFORM:
        for (int i=-emf->B.gc[0]; i<emf->B.nx+emf->B.gc[1]; i++) {
            Bx[ i ] = init_fld -> B_0.x;
            By[ i ] = init_fld -> B_0.y;
            Bz[ i ] = init_fld -> B_0.z;
        }
        fftr_r2c( emf->fft_forward, Bx, emf->fB.x );
        fftr_r2c( emf->fft_forward, By, emf->fB.y );
        fftr_r2c( emf->fft_forward, Bz, emf->fB.z );

        break;

    case EMF_FLD_TYPE_CUSTOM:
        for (int i=-emf->B.gc[0]; i<emf->B.nx+emf->B.gc[1]; i++) {
            float3 init_B = (init_fld->B_custom)
                (i,emf->dx, init_fld->B_custom_data);
            Bx[ i ] = init_B.x;
            By[ i ] = init_B.y;
            Bz[ i ] = init_B.z;
        }
        fftr_r2c( emf->fft_forward, Bx, emf->fB.x );
        fftr_r2c( emf->fft_forward, By, emf->fB.y );
        fftr_r2c( emf->fft_forward, Bz, emf->fB.z );
        break;
    }    

}

