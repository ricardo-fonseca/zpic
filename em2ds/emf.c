/**
 * @file emf.c
 * @author Ricardo Fonseca
 * @brief EM fields
 * @version 0.2
 * @date 2022-02-04
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
 * @brief Initializes the EM field objecnt
 * 
 * @param emf			EM fields object
 * @param nx			Number of cells
 * @param box			Physical box size
 * @param dt			Simulation time step
 */
void emf_new( t_emf *emf, const int nx[], float box[], const float dt )
{

	if ( nx[0] % 2 || nx[1] % 2 ) {
		fprintf(stderr,"(*error*) Only even grid sizes are supported.");
		exit(-1);
	}

	// Number of guard cells for linear interpolation
	const unsigned int gc[2][2] = {{0,1},
		                           {0,1}};

	// The FFT result will be transposed
	const unsigned int fnx[2] = { nx[1], nx[0]/2+1 };

	// Initialize grids
	float3_grid2d_init( &emf->E, (unsigned int *) nx, gc );
	float3_grid2d_init( &emf->B, (unsigned int *) nx, gc );

	cfloat3_grid2d_init( &emf->fEl, fnx, NULL );
	cfloat3_grid2d_init( &emf->fEt, fnx, NULL );
	cfloat3_grid2d_init( &emf->fB,  fnx, NULL );

	// zero fields
	float3_grid2d_zero( &emf->E );
	float3_grid2d_zero( &emf->B );

	cfloat3_grid2d_zero( &emf->fEl );
	cfloat3_grid2d_zero( &emf->fEt );
	cfloat3_grid2d_zero( &emf->fB  );

	// Initializ FFT transforms
	fftr2d_init_cfg( &emf -> fft_forward,  nx[0], nx[1], emf -> E.nrow, FFT_FORWARD );
	fftr2d_init_cfg( &emf -> fft_backward, nx[0], nx[1], emf -> E.nrow, FFT_BACKWARD );

	// Set cell sizes and box limits
	for(int i = 0; i<2; i++){
		emf -> box[i] = box[i];
		emf -> dx[i] = box[i] / nx[i];
	}

    // store time step values
    emf -> dt = dt;

	// Reset iteration number
	emf -> iter = 0;

	// Set default solver type to Pseudo-spectral Analytic Time Domain
	emf -> solver_type = EMF_SOLVER_PSATD;
//	emf -> solver_type = EMF_SOLVER_PSTD;

}

/**
 * @brief Frees dynamic memory from EM fields.
 * 
 * @param emf 	EM fields
 */
void emf_delete( t_emf *emf )
{
	float3_grid2d_cleanup( &emf->E );
	float3_grid2d_cleanup( &emf->B );

	cfloat3_grid2d_cleanup( &emf->fEl );
	cfloat3_grid2d_cleanup( &emf->fEt );
	cfloat3_grid2d_cleanup( &emf->fB );

	fftr2d_cleanup_cfg( &emf -> fft_backward );
	fftr2d_cleanup_cfg( &emf -> fft_forward );
}

/*********************************************************************************************

 Laser Pulses

*********************************************************************************************/

/**
 * @brief Gaussian laser pulse phase and transverse envelope
 * 
 * @param laser 	Laser pulse parameters
 * @param z_ 		Longitudinal position
 * @param r 		Radial position (distance from propagation axis)
 * @return 			Gaussian phase
 */
float gauss_phase( const t_emf_laser* const laser, const float z_, const float r )
{
	float z = z_ - laser -> focus;
	float z0   = laser -> omega0 * ( laser->W0 * laser->W0 ) / 2;
	float rho2 = r*r;
	float curv = rho2 * z / (z0*z0 + z*z);
	float rWl2 = (z0*z0)/(z0*z0 + z*z);
	float gouy_shift = atan2( z, z0 );

	return sqrt( sqrt(rWl2) ) *
	exp( - rho2 * rWl2/(laser->W0 * laser->W0) ) *
	cos( laser -> omega0*( z + curv ) - gouy_shift );
}

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
 * @brief Ensures div E = div B = 0 for the whole grid by modifying the 
 * 		x component of E and B
 * 
 * This is used to calculate the longitudinal fields associated with a 
 * laser pulse
 * 
 * @param fld 		E or B fields
 * @param dk		k space cell size
 */
void div_corr_x( t_cfloat3_grid2d *fld, const float dk[] )
{
    float complex * const restrict fldx = fld->x;
    float complex * const restrict fldy = fld->y;

	const int fnrow = fld -> nrow;

	// Calculates the x component of the field by enforcing div(fld) = 0

	// no work is required for kx = 0
	for ( int i = 1; i < fld->nx[1]; i++) {
		float kx = i * dk[0];
		for( int j = 0; j < fld->nx[0]; j++) {
			float ky = ((j <= fld->nx[0]/2) ? j : (j - (int) fld->nx[0]) ) * dk[1];

			fldx[ i * fnrow + j ] = - ky * fldy[ i * fnrow + j ] / kx;
		}
	}
}

/**
 * @brief Add laser pulse to simulation.
 * 
 * Laser pulses are superimposed on top of existing E and B fields. 
 * Currently only 1 gaussian pulse can be added
 * 
 * @param emf 		EM fields
 * @param laser 	Laser pulse parameters
 */
void emf_add_laser( t_emf* const emf, const t_emf_laser* const laser )
{
	int i, j, nrow;

	float z_center, r_center, z, r;
	float amp, lenv, k;
	float dx, dy;
	float cos_pol, sin_pol;

	float* Ey = (float *) malloc( emf->E.nx[0] * emf->E.nx[1] * sizeof(float));
	float* Ez = (float *) malloc( emf->E.nx[0] * emf->E.nx[1] * sizeof(float));
	float* By = (float *) malloc( emf->E.nx[0] * emf->E.nx[1] * sizeof(float));
	float* Bz = (float *) malloc( emf->E.nx[0] * emf->E.nx[1] * sizeof(float));

	nrow = emf -> E.nx[0];
	dx = emf -> dx[0];
	dy = emf -> dx[1];

	z_center = laser->start - laser->fwhm/2;
	r_center = laser->axis;
	amp = laser->omega0 * laser->a0;

	cos_pol = cos( laser -> polarization );
	sin_pol = sin( laser -> polarization );

	switch (laser->type) {
		case PLANE:
			k = laser -> omega0;

			for (i = 0; i < emf->E.nx[0]; i++) {
				z = i * dx - z_center;

				lenv   = amp*lon_env( laser, z );

				for (j = 0; j < emf->E.nx[1]; j++) {
					// E[i + j*nrow].x += 0.0
					Ey[i + j*nrow] = +lenv * cos( k * z ) * cos_pol;
					Ez[i + j*nrow] = +lenv * cos( k * z ) * sin_pol;

					// E[i + j*nrow].x += 0.0
					By[i + j*nrow] = -lenv * cos( k * z ) * sin_pol;
					Bz[i + j*nrow] = +lenv * cos( k * z ) * cos_pol;

				}
			}
			break;

		case GAUSSIAN:

			for (i = 0; i < emf->E.nx[0]; i++) {
				z = i * dx - z_center;

				lenv   = amp*lon_env( laser, z );

				for (j = 0; j < emf->E.nx[1]; j++) {
					r = j * dy - r_center;

					// E[i + j*nrow].x += 0.0
					Ey[i + j*nrow] = +lenv * gauss_phase( laser, z, r ) * cos_pol;
					Ez[i + j*nrow] = +lenv * gauss_phase( laser, z, r ) * sin_pol;

					// B[i + j*nrow].x += 0.0
					By[i + j*nrow] = -lenv * gauss_phase( laser, z, r ) * sin_pol;
					Bz[i + j*nrow] = +lenv * gauss_phase( laser, z, r ) * cos_pol;

				}
			}

			break;
		default:
			break;
	}

	// Buffer to hold Fourier transform of laser field
	t_cfloat3_grid2d buffer;
	cfloat3_grid2d_init( &buffer, (unsigned int *) emf->fEt.nx, NULL);
	cfloat3_grid2d_zero( &buffer );

	t_fftr2d_cfg fft_forward;
	fftr2d_init_cfg( &fft_forward, emf->E.nx[0], emf->E.nx[1],
	                 0, FFT_FORWARD );

	const int nrowfEt = emf->fEt.nrow;
	const int nrowbuf = emf->fEt.nx[0];

	const float dk[2] = { fft_dk( emf->E.nx[0], emf->dx[0] ),
		                  fft_dk( emf->E.nx[1], emf->dx[1] ) };

	const int kx1 = (emf->fEt.nx[1]-1)/2;
	const int ky1 = emf->fEt.nx[0]/4;
	const int ky2 = emf->fEt.nx[0] - ky1;

	// Transform transverse components of laser E-field
	fftr2d_r2c( &fft_forward, Ey, buffer.y );
	fftr2d_r2c( &fft_forward, Ez, buffer.z );

	if ( laser->type != PLANE ) {
		// Solve for Ex using div E = 0
		div_corr_x( &buffer, dk );
	}

	// Add to transverse E field
	for(i = 0; i < kx1; i++) {
		for( j = 0; j < ky1; j++) {
			emf->fEt.x[ i * nrowfEt + j ] += buffer.x[ i * nrowbuf + j ];
			emf->fEt.y[ i * nrowfEt + j ] += buffer.y[ i * nrowbuf + j ];
			emf->fEt.z[ i * nrowfEt + j ] += buffer.z[ i * nrowbuf + j ];
		}
		for( j = ky2; j < emf->fEt.nx[0]; j++) {
			emf->fEt.x[ i * nrowfEt + j ] += buffer.x[ i * nrowbuf + j ];
			emf->fEt.y[ i * nrowfEt + j ] += buffer.y[ i * nrowbuf + j ];
			emf->fEt.z[ i * nrowfEt + j ] += buffer.z[ i * nrowbuf + j ];
		}
	}

	// Transform transverse components of laser B-field
	fftr2d_r2c( &fft_forward, By, buffer.y );
	fftr2d_r2c( &fft_forward, Bz, buffer.z );

	if ( laser->type != PLANE ) {
		// Solve for Ex using div E = 0
		div_corr_x( &buffer, dk );
	}

	for(i = 0; i < kx1; i++) {
		for( j = 0; j < ky1; j++) {
			emf->fB.x[ i * nrowfEt + j ] += buffer.x[ i * nrowbuf + j ];
			emf->fB.y[ i * nrowfEt + j ] += buffer.y[ i * nrowbuf + j ];
			emf->fB.z[ i * nrowfEt + j ] += buffer.z[ i * nrowbuf + j ];
		}
		for( j = ky2; j < emf->fEt.nx[0]; j++) {
			emf->fB.x[ i * nrowfEt + j ] += buffer.x[ i * nrowbuf + j ];
			emf->fB.y[ i * nrowfEt + j ] += buffer.y[ i * nrowbuf + j ];
			emf->fB.z[ i * nrowfEt + j ] += buffer.z[ i * nrowbuf + j ];
		}
	}

	// Cleanup fft configuration
	fftr2d_cleanup_cfg( &fft_forward );

	// Cleanup temporary values
	cfloat3_grid2d_cleanup( &buffer );

	free( Bz );
	free( By );
	free( Ez );
	free( Ey );


	// Transform to real fields for diagnostics
   	fftr2d_c2r( &emf -> fft_backward, emf -> fEt.x, emf ->E.x );
   	fftr2d_c2r( &emf -> fft_backward, emf -> fEt.y, emf ->E.y );
   	fftr2d_c2r( &emf -> fft_backward, emf -> fEt.z, emf ->E.z );

   	fftr2d_c2r( &emf -> fft_backward, emf -> fB.x, emf ->B.x );
   	fftr2d_c2r( &emf -> fft_backward, emf -> fB.y, emf ->B.y );
   	fftr2d_c2r( &emf -> fft_backward, emf -> fB.z, emf ->B.z );
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
	const t_float3_grid2d * vfld;
	switch (field) {
		case EFLD:
			vfld = &emf->E;
            snprintf(vfname,16,"E%1d",fc);
            snprintf(vflabel,16,"E_%c",comp[fc]);
			break;
		case BFLD:
			vfld = &emf->B;
            snprintf(vfname,16,"B%1d",fc);
            snprintf(vflabel,16,"B_%c",comp[fc]);
			break;
		default:
			fprintf(stderr,"(*error*) Invalid field type selected, returning\n");
			return;
	}

	// Pack the information
	float * restrict const buf = malloc( vfld->nx[0]*vfld->nx[1]*sizeof(float) );
    float * restrict p = buf;
    float * f;

	switch (fc) {
		case 0:
			f = vfld->x;
			break;
		case 1:
			f = vfld->y;
			break;
		case 2:
			f = vfld->z;
			break;
		default:
			fprintf(stderr,"(*error*) Invalid field component selected, returning\n");
			return;
	}

	for( int j = 0; j < vfld->nx[1]; j++) {
		for ( int i = 0; i < vfld->nx[0]; i++ ) {
			p[i] = f[i];
		}
		p += vfld->nx[0];
		f += vfld->nrow;
	}


    t_zdf_grid_axis axis[2];
    axis[0] = (t_zdf_grid_axis) {
    	.min = 0.0,
    	.max = emf->box[0],
		.name = "x",
    	.label = "x",
    	.units = "c/\\omega_p"
    };

    axis[1] = (t_zdf_grid_axis) {
    	.min = 0.0,
    	.max = emf->box[1],
		.name = "y",
    	.label = "y",
    	.units = "c/\\omega_p"
    };

    t_zdf_grid_info info = {
    	.ndims = 2,
		.name = vfname,
    	.label = vflabel,
    	.units = "m_e c \\omega_p e^{-1}",
    	.axis = axis
    };

    info.count[0] = vfld->nx[0];
    info.count[1] = vfld->nx[1];

    t_zdf_iteration iter = {
		.name = "ITERATION",
    	.n = emf->iter,
    	.t = emf -> iter * emf -> dt,
    	.time_units = "1/\\omega_p"
    };

	zdf_save_grid( (float *) buf, zdf_float32, &info, &iter, "EMF" );

	// free local data
	free( buf );

}


/*********************************************************************************************

 Pseudo-spectral Time Domain Field solver

 *********************************************************************************************/

/**
 * @brief Advance the magnetic field using a PSTD algorithm
 * 
 * @param emf   EM fields
 * @param dt    Time step
 */
void advance_fB( t_emf *emf, const float dt )
{
    float complex * const restrict fBx = emf -> fB.x;
    float complex * const restrict fBy = emf -> fB.y;
    float complex * const restrict fBz = emf -> fB.z;

    float complex * const restrict fEtx = emf -> fEt.x;
    float complex * const restrict fEty = emf -> fEt.y;
    float complex * const restrict fEtz = emf -> fEt.z;

    const float dkx = fft_dk( emf->E.nx[0], emf->dx[0] );
    const float dky = fft_dk( emf->E.nx[1], emf->dx[1] );

	const int fnrow = emf->fB.nrow;

	// Canonical implementation
	for (int i = 0; i < emf -> fB.nx[1]; i++) {
		float kx = i * dkx;
		for (int j = 0; j < emf -> fB.nx[0]; j++) {
			float ky = ((j <= emf -> fB.nx[0]/2) ? j : ( j - (int) emf -> fB.nx[0] ) ) * dky;
			const unsigned int idx = i * fnrow + j;

			fBx[idx] +=  - dt * I * ky * fEtz[idx];
			fBy[idx] +=  + dt * I * kx * fEtz[idx];
			fBz[idx] +=  + dt * I * ( -kx * fEty[idx] + ky * fEtx[idx] );
		}
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

	float complex * const restrict fEtx = emf -> fEt.x;
	float complex * const restrict fEty = emf -> fEt.y;
	float complex * const restrict fEtz = emf -> fEt.z;

	float complex * const restrict fBx = emf -> fB.x;
	float complex * const restrict fBy = emf -> fB.y;
	float complex * const restrict fBz = emf -> fB.z;

	float complex * const restrict fJx = current -> fJ.x;
	float complex * const restrict fJy = current -> fJ.y;
	float complex * const restrict fJz = current -> fJ.z;

	const float dkx = fft_dk( emf->E.nx[0], emf->dx[0] );
	const float dky = fft_dk( emf->E.nx[1], emf->dx[1] );

	const int fnrow = emf->fEt.nrow;

	// kx = 0, ky = 0
	/*
	Note that we are zeroing the x and y components for k = 0
	 */

	{
		// const float complex fJtx = 0;
		// const float complex fJty = 0;
		const float complex fJtz = fJz[0];

		// fEtx[0] -= dt * fJtx;
		// fEty[0] -= dt * fJty;
		fEtz[0] -= dt * fJtz;
	}

	// kx = 0, ky != 0
	for (int j = 1; j < emf -> fEt.nx[0]; j++) {
		const float ky = ((j <= emf -> fEt.nx[0]/2) ? j : (j - (int) emf -> fEt.nx[0]) ) * dky;

		// Get transverse current
		const float complex fJtx = fJx[j];
		const float complex fJty = 0;
		const float complex fJtz = fJz[j];

		// Advance field
		fEtx[j] += dt * ( +I * ky * fBz[j] - fJtx );
		fEty[j] += dt * (                  - fJty );
		fEtz[j] += dt * ( -I * ky * fBx[j] - fJtz );
	}


	// kx = 0, ky != 0
	for (int i = 1; i < emf -> fEt.nx[1]; i++) {
		const float kx = i * dkx;
		for (int j = 0; j < emf -> fEt.nx[0]; j++) {
			const float ky = ((j <= emf -> fEt.nx[0]/2) ? j : (j - (int) emf -> fEt.nx[0]) ) * dky;
			const unsigned int idx = i * fnrow + j;

			// Get transverse current
			const float complex kdJ_k2 = (kx * fJx[idx] + ky * fJy[idx]) / (kx*kx + ky*ky);
			const float complex fJtx = fJx[idx] - kx * kdJ_k2;
			const float complex fJty = fJy[idx] - ky * kdJ_k2;
			const float complex fJtz = fJz[idx];

			// Advance field
			fEtx[idx] += dt * ( +I *   ky * fBz[idx]                   - fJtx );
			fEty[idx] += dt * ( -I *   kx * fBz[idx]                   - fJty );
			fEtz[idx] += dt * ( +I * ( kx * fBy[idx] - ky * fBx[idx] ) - fJtz );

		}
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

	const float dkx = fft_dk( emf->E.nx[0], emf->dx[0] );
	const float dky = fft_dk( emf->E.nx[1], emf->dx[1] );

	const int fnrow = emf->fEl.nrow;

	// Calculates the longitudinal component of E from the
	// charge density:
	// $E_L(\mathbf{k}) = - i \frac{ \mathbf{k} \rho(\mathbf{k})}{k^2}$

	// Take care of special case k = 0
	fElx[0] = 0;
	fEly[0] = 0;
	fElz[0] = 0;

	// Take care of special cases kx = 0
	for (int j = 1; j < emf -> fEl.nx[0]; j++) {
		const float ky = ((j <= emf -> fEl.nx[0]/2) ? j : (j - (int) emf -> fEl.nx[0]) ) * dky;

		fElx[j] = 0;
		fEly[j] = -I * frho[j] / ky;
		fElz[j] = 0 ;
	}

	// Solve remaining cases
	for (int i = 1; i < emf -> fEl.nx[1]; i++) {
		float kx = i * dkx;
		for (int j = 0; j < emf -> fEl.nx[0]; j++) {
			const float ky = ((j <= emf -> fEl.nx[0]/2) ? j : (j - (int) emf -> fEl.nx[0] )) * dky;
			const float k2 = kx*kx + ky*ky;

			const unsigned int idx = i * fnrow + j;

			fElx[idx] = -I * kx * frho[idx] / k2;
			fEly[idx] = -I * ky * frho[idx] / k2;
			fElz[idx] = 0 ;
		}
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

    float complex * const restrict fEtx = emf -> fEt.x;
    float complex * const restrict fEty = emf -> fEt.y;
    float complex * const restrict fEtz = emf -> fEt.z;

    float complex * const restrict fBx = emf -> fB.x;
    float complex * const restrict fBy = emf -> fB.y;
    float complex * const restrict fBz = emf -> fB.z;

    // Full current
    float complex * const restrict fJx = current -> fJ.x;
    float complex * const restrict fJy = current -> fJ.y;
    float complex * const restrict fJz = current -> fJ.z;

	const float dkx = fft_dk( emf->E.nx[0], emf->dx[0] );
	const float dky = fft_dk( emf->E.nx[1], emf->dx[1] );

	// This is the same value for fEt, fB and fJ
	const int fnrow = emf->fEt.nrow;

	// kx = ky = 0
	{
		/*
		 Note that we are setting fJtx = fJty = 0 for k = 0, meaning that we are eliminating any
		 net transverse current along x and y.
		 */

		// const float complex fJtx = 0;
		// const float complex fJty = 0;
		const float complex fJtz = fJz[0];

		//fEtx[0] -= dt * fJtx;		// Not needed because fJtx = 0
		//fEty[0] -= dt * fJty;		// Not needed because fJty = 0
		fEtz[0] -= dt * fJtz;

		// fBx[0]  += 0;
		// fBy[0]  += 0;
		// fBz[0]  += 0;
	}

	// kx = 0, ky != 0
	for (int j = 1; j < emf -> fEt.nx[0]; j++) {
		const float ky = ((j <= emf -> fEt.nx[0]/2) ? j : (j - (int) emf -> fEt.nx[0]) ) * dky;
		const float k2 = ky*ky;

		// Calculate transverse current
		const float complex fJtx = fJx[j];
		const float complex fJty = 0;
		const float complex fJtz = fJz[j];

		const float C   = cosf( ky * dt );
		const float S_k = sinf( ky * dt ) / ky;
		const float complex I1mC_k2 = I * (1.0f - C) / k2;

		// PSATD Field advance equations
		float complex Ex = fEtx[j];
		float complex Ey = fEty[j];
		float complex Ez = fEtz[j];

		float complex Bx = fBx[j];
		float complex By = fBy[j];
		float complex Bz = fBz[j];

		Ex = C * Ex + S_k * (  I * ky * fBz[j] - fJtx );
		Ey = C * Ey + S_k * (                  - fJty );
		Ez = C * Ez + S_k * ( -I * ky * fBx[j] - fJtz );

		Bx = C * Bx - S_k * I * ky * fEtz[j] + I1mC_k2 * ky * fJz[j];
		By = C * By;
		Bz = C * Bz + S_k * I * ky * fEtx[j] - I1mC_k2 * ky * fJx[j];

		fEtx[j] = Ex;
		fEty[j] = Ey;
		fEtz[j] = Ez;

		fBx[j]  = Bx;
		fBy[j]  = By;
		fBz[j]  = Bz;

	}


	// kx > 0
	for (int i = 1; i < emf -> fEt.nx[1]; i++) {
		const float kx  = i * dkx;

		for (int j = 0; j < emf -> fEt.nx[0]; j++) {
			const float ky = ((j <= emf -> fEt.nx[0]/2) ? j : (j - (int) emf -> fEt.nx[0]) ) * dky;

			const float k2 = kx*kx + ky*ky;
			const float k = sqrtf(k2);

			const unsigned int idx = i * fnrow + j;

			// Calculate transverse current
			const float complex kdJ_k2 = (kx * fJx[idx] + ky * fJy[idx])/k2;
			const float complex fJtx = fJx[idx] - kx * kdJ_k2;
			const float complex fJty = fJy[idx] - ky * kdJ_k2;
			const float complex fJtz = fJz[idx];

			// PSATD Field advance equations
			const float C   = cosf( k * dt );
			const float S_k = sinf( k * dt ) / k;
			const float complex I1mC_k2 = I * (1.0f - C) / k2;

			float complex Ex = fEtx[idx];
			float complex Ey = fEty[idx];
			float complex Ez = fEtz[idx];

			float complex Bx = fBx[idx];
			float complex By = fBy[idx];
			float complex Bz = fBz[idx];

			Ex = C * Ex + S_k * ( I * (  ky *  fBz[idx]                  ) - fJtx );
			Ey = C * Ey + S_k * ( I * ( -kx *  fBz[idx]                  ) - fJty );
			Ez = C * Ez + S_k * ( I * (  kx *  fBy[idx] - ky *  fBx[idx] ) - fJtz );

			Bx = C * Bx - S_k * ( I * (  ky * fEtz[idx]                  ) ) + I1mC_k2 * (  ky * fJz[idx]                 );
			By = C * By - S_k * ( I * ( -kx * fEtz[idx]                  ) ) + I1mC_k2 * ( -kx * fJz[idx]                 );
			Bz = C * Bz - S_k * ( I * (  kx * fEty[idx] - ky * fEtx[idx] ) ) + I1mC_k2 * (  kx * fJy[idx] - ky * fJx[idx] );

			fEtx[idx] = Ex;
			fEty[idx] = Ey;
			fEtz[idx] = Ez;

			fBx[idx]  = Bx;
			fBy[idx]  = By;
			fBz[idx]  = Bz;

		}
	}
}


/**
 * @brief Updates real electric and magnetic field from Fourier transforms
 * 
 * It will also update the guard cell values for field interpolation
 * 
 * @param emf   EM fields
 */
void emf_update( t_emf *emf )
{
	// Update E field

	// Add transverse and longitudinal components
	cfloat3_grid2d_add( &emf -> fEl, &emf -> fEt );

   	// Transform to real fields
   	fftr2d_c2r( &emf -> fft_backward, emf -> fEl.x, emf -> E.x );
   	fftr2d_c2r( &emf -> fft_backward, emf -> fEl.y, emf -> E.y );
   	fftr2d_c2r( &emf -> fft_backward, emf -> fEl.z, emf -> E.z );

   	// Update B field
   	fftr2d_c2r( &emf -> fft_backward, emf -> fB.x, emf -> B.x );
   	fftr2d_c2r( &emf -> fft_backward, emf -> fB.y, emf -> B.y );
   	fftr2d_c2r( &emf -> fft_backward, emf -> fB.z, emf -> B.z );

	// Update guard cells
	const int nrow = emf->E.nrow;

	float* restrict const Ex = emf -> E.x;
	float* restrict const Ey = emf -> E.y;
	float* restrict const Ez = emf -> E.z;

	float* restrict const Bx = emf -> B.x;
	float* restrict const By = emf -> B.y;
	float* restrict const Bz = emf -> B.z;

	int const nx0  = emf -> E.nx[0];
	int const nx1  = emf -> E.nx[1];

	int const gc00 = emf -> E.gc[0][0];
	int const gc01 = emf -> E.gc[0][1];
	int const gc10 = emf -> E.gc[1][0];
	int const gc11 = emf -> E.gc[1][1];

	// x
	for (int j = -gc10; j < nx1 + gc11; j++) {

		// lower
		for (int i = -gc00; i < 0; i++) {
			Ex[ i + j*nrow ] = Ex[ nx0 + i + j*nrow ];
			Ey[ i + j*nrow ] = Ey[ nx0 + i + j*nrow ];
			Ez[ i + j*nrow ] = Ez[ nx0 + i + j*nrow ];

			Bx[ i + j*nrow ] = Bx[ nx0 + i + j*nrow ];
			By[ i + j*nrow ] = By[ nx0 + i + j*nrow ];
			Bz[ i + j*nrow ] = Bz[ nx0 + i + j*nrow ];
		}

		// upper
		for (int i = 0; i < gc01; i++) {
			Ex[ nx0 + i + j*nrow ] = Ex[ i + j*nrow ];
			Ey[ nx0 + i + j*nrow ] = Ey[ i + j*nrow ];
			Ez[ nx0 + i + j*nrow ] = Ez[ i + j*nrow ];

			Bx[ nx0 + i + j*nrow ] = Bx[ i + j*nrow ];
			By[ nx0 + i + j*nrow ] = By[ i + j*nrow ];
			Bz[ nx0 + i + j*nrow ] = Bz[ i + j*nrow ];
		}

	}


	// y
	for (int i = -gc00; i < nx0 + gc01; i++) {

		// lower
		for (int j=-gc10; j<0; j++) {
			Ex[ i + j*nrow ] = Ex[ i + (nx1+j)*nrow ];
			Ey[ i + j*nrow ] = Ey[ i + (nx1+j)*nrow ];
			Ez[ i + j*nrow ] = Ez[ i + (nx1+j)*nrow ];

			Bx[ i + j*nrow ] = Bx[ i + (nx1+j)*nrow ];
			By[ i + j*nrow ] = By[ i + (nx1+j)*nrow ];
			Bz[ i + j*nrow ] = Bz[ i + (nx1+j)*nrow ];
		}

		// upper
		for (int j=0; j<gc11; j++) {
			Ex[ i + (nx1+j)*nrow ] = Ex[ i + j*nrow ];
			Ey[ i + (nx1+j)*nrow ] = Ey[ i + j*nrow ];
			Ez[ i + (nx1+j)*nrow ] = Ez[ i + j*nrow ];

			Bx[ i + (nx1+j)*nrow ] = Bx[ i + j*nrow ];
			By[ i + (nx1+j)*nrow ] = By[ i + j*nrow ];
			Bz[ i + (nx1+j)*nrow ] = Bz[ i + j*nrow ];
		}
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

	float* restrict const Ex = emf -> E.x;
	float* restrict const Ey = emf -> E.y;
	float* restrict const Ez = emf -> E.z;

	float* restrict const Bx = emf -> B.x;
	float* restrict const By = emf -> B.y;
	float* restrict const Bz = emf -> B.z;

    const int nrow = emf -> E.nrow;

	for( int i = 0; i<6; i++) energy[i] = 0;

	for( int j = 0; j < emf -> E.nx[1]; j ++ ) {
		for( int i = 0; i < emf -> E.nx[0]; i ++ ) {
			energy[0] += Ex[i + j*nrow] * Ex[i + j*nrow];
			energy[1] += Ey[i + j*nrow] * Ey[i + j*nrow];
			energy[2] += Ez[i + j*nrow] * Ez[i + j*nrow];
			energy[3] += Bx[i + j*nrow] * Bx[i + j*nrow];
			energy[4] += By[i + j*nrow] * By[i + j*nrow];
			energy[5] += Bz[i + j*nrow] * Bz[i + j*nrow];
		}
	}

	for( int i = 0; i<6; i++) energy[i] *= 0.5 * emf -> dx[0] * emf -> dx[1];
}
