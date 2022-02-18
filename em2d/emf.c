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

void emf_move_window( t_emf *emf );
void emf_update_part_fld( t_emf *emf );
void emf_update_gc( t_emf *emf );

/// Time spent advancing the EM fields
static double _emf_time = 0.0;

/**
 * @brief Time spent advancing the EM fields
 * 
 * @return 		Time spent in seconds
 */
double emf_time( void )
{
	return _emf_time;
}

/*********************************************************************************************
 
 Constructor / Destructor
 
 *********************************************************************************************/

/**
 * @brief Initalized EM fields object
 * 
 * Fields are initialized with 0 values, if you require other initial
 * values use the `init_fld()` function.
 * 
 * @param emf 	EM fields
 * @param nx 	Number of grid cells [x,y]
 * @param box 	Physical box size [x,y]
 * @param dt 	Simulation time step
 */
void emf_new( t_emf *emf, int nx[], float box[], const float dt )
{   
	// Number of guard cells for linear interpolation
	int gc[2][2] = {{1,2},
		            {1,2}}; 
	
	// Allocate global arrays
	size_t size;
	
	size = (gc[0][0] + nx[0] + gc[0][1]) * (gc[1][0] + nx[1] + gc[1][1]) * sizeof( float3 );
	
	emf->E_buf = malloc( size );
	emf->B_buf = malloc( size );
	
	assert( emf->E_buf && emf->B_buf );
	
	// zero fields
	memset( emf->E_buf, 0, size );
	memset( emf->B_buf, 0, size );
	
	// store nx and gc values
	for(int i = 0; i<2; i++){
		emf->nx[i] = nx[i];
		emf->gc[i][0] = gc[i][0];
		emf->gc[i][1] = gc[i][1];
	}
	emf -> nrow = gc[0][0] + nx[0] + gc[0][1];
	
    // store time step values
    emf -> dt = dt;

	// Make E and B point to cell [0][0]
	emf->E = emf->E_buf + gc[0][0] + gc[1][0] * emf->nrow;
	emf->B = emf->B_buf + gc[0][0] + gc[1][0] * emf->nrow;
	
	// Set cell sizes and box limits
	for(int i = 0; i<2; i++){
		emf -> box[i] = box[i];
		emf -> dx[i] = box[i] / nx[i];
	}

	// Set time step
	emf -> dt = dt;

	// Reset iteration number
	emf -> iter = 0;

	// Reset moving window information
	emf -> moving_window = 0;
	emf -> n_move = 0;

    // Disable external fields by default
	emf -> ext_fld.E_type = EMF_FLD_TYPE_NONE;
	emf -> ext_fld.B_type = EMF_FLD_TYPE_NONE;
	emf -> E_part = emf->E;
	emf -> B_part = emf->B;

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
	free( emf->E_buf );
	free( emf->B_buf );
	
	emf->E_buf = NULL;
	emf->B_buf = NULL;

	if ( emf -> ext_fld.E_type > EMF_FLD_TYPE_NONE ) {
		free( emf -> ext_fld.E_part_buf );
	}

	if ( emf -> ext_fld.B_type > EMF_FLD_TYPE_NONE ) {
		free( emf -> ext_fld.B_part_buf );
	}

	emf->E_part = NULL;
	emf->B_part = NULL;
	
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
 * @return 			Envelope value
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
 * @param emf 		EM fields
 */
void div_corr_x( t_emf *emf )
{
	float3* restrict E = emf -> E;
	float3* restrict B = emf -> B;
	const int nrow = emf -> nrow;
	const double dx_dy = emf -> dx[0]/emf -> dx[1];
	
	for (int j=0; j<emf->nx[1]; j++) {
		double ex = 0.0;
		double bx = 0.0;
		for (int i = emf->nx[0]-1; i>=0; i--) {
			ex += dx_dy * (E[i+1 + j*nrow].y - E[i+1 + (j-1)*nrow ].y);
			E[i+j*nrow].x = ex;
			
			bx += dx_dy * (B[i + (j+1)*nrow].y - B[i + j*nrow ].y);
			B[i+j*nrow].x = bx;
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
void emf_add_laser( t_emf* const emf,  t_emf_laser*  laser )
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
	int i, j, nrow;
	
	float r_center, z, z_2, r, r_2;
	float amp, lenv, lenv_2, k;
	float dx, dy;
	float cos_pol, sin_pol;
	
	float3* restrict E = emf -> E;
	float3* restrict B = emf -> B;

	nrow = emf -> nrow;
	dx = emf -> dx[0];
	dy = emf -> dx[1];
	
	r_center = laser->axis;
	amp = laser->omega0 * laser->a0;
	
	cos_pol = cos( laser -> polarization );
	sin_pol = sin( laser -> polarization );
		
	switch (laser->type) {
		case PLANE:
			k = laser -> omega0;

			for (i = 0; i < emf->nx[0]; i++) {
				z = i * dx;
				z_2 = z + dx/2;
				
				lenv   = amp*lon_env( laser, z );
				lenv_2 = amp*lon_env( laser, z_2 );
				
				for (j = 0; j < emf->nx[1]; j++) {
					// E[i + j*nrow].x += 0.0
					E[i + j*nrow].y += +lenv * cos( k * z ) * cos_pol;
					E[i + j*nrow].z += +lenv * cos( k * z ) * sin_pol;

					// E[i + j*nrow].x += 0.0
					B[i + j*nrow].y += -lenv_2 * cos( k * z_2 ) * sin_pol;
					B[i + j*nrow].z += +lenv_2 * cos( k * z_2 ) * cos_pol;
					
				}
			}
			break;

		case GAUSSIAN:

			for (i = 0; i < emf->nx[0]; i++) {
				z = i * dx;
				z_2 = z + dx/2;
				
				lenv   = amp*lon_env( laser, z );
				lenv_2 = amp*lon_env( laser, z_2 );
				
				for (j = 0; j < emf->nx[1]; j++) {
					r = j * dy - r_center;
					r_2 = r + dy/2;
										
					// E[i + j*nrow].x += 0.0
					E[i + j*nrow].y += +lenv * gauss_phase( laser, z  , r_2 ) * cos_pol;
					E[i + j*nrow].z += +lenv * gauss_phase( laser, z  , r   ) * sin_pol;
					
					// B[i + j*nrow].x += 0.0
					B[i + j*nrow].y += -lenv_2 * gauss_phase( laser, z_2, r   ) * sin_pol;
					B[i + j*nrow].z += +lenv_2 * gauss_phase( laser, z_2, r_2 ) * cos_pol;
					
				}
			}
			div_corr_x(emf);
			
			break;
		default:
			break;
	}
    
	// Set guard cell values
	emf_update_gc( emf );
	
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
	float3 * restrict f;
	switch (field) {
		case EFLD:
			f = emf->E;
            snprintf(vfname,16,"E%1d",fc);
            snprintf(vflabel,16,"E_%c",comp[fc]);
			break;
		case BFLD:
			f = emf->B;
            snprintf(vfname,16,"B%1d",fc);
            snprintf(vflabel,16,"B_%c",comp[fc]);
			break;
		case EPART:
			f = emf->E_part;
            snprintf(vfname,16,"E%1d-part",fc);
            snprintf(vflabel,16,"E_{%cp}",comp[fc]);
			break;
		case BPART:
			f = emf->B_part;
            snprintf(vfname,16,"B%1d-part",fc);
            snprintf(vflabel,16,"B_{%cp}",comp[fc]);
			break;
		default:
			fprintf(stderr, "Invalid field type selected, returning\n");
			return;
	}

	// Pack the information
	float * restrict const buf = malloc( emf->nx[0]*emf->nx[1]*sizeof(float) );
    float * restrict p = buf;
	switch (fc) {
		case 0:
			for( int j = 0; j < emf->nx[1]; j++) {
				for ( int i = 0; i < emf->nx[0]; i++ ) {
					p[i] = f[i].x;
				}
				p += emf->nx[0];
				f += emf->nrow;
			}
			break;
		case 1:
			for( int j = 0; j < emf->nx[1]; j++) {
				for ( int i = 0; i < emf->nx[0]; i++ ) {
					p[i] = f[i].y;
				}
				p += emf->nx[0];
				f += emf->nrow;
			}
			break;
		case 2:
			for( int j = 0; j < emf->nx[1]; j++) {
				for ( int i = 0; i < emf->nx[0]; i++ ) {
					p[i] = f[i].z;
				}
				p += emf->nx[0];
				f += emf->nrow;
			}
			break;
		default:
			fprintf(stderr, "Invalid field component selected, returning\n");
			return;
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

    info.count[0] = emf->nx[0];
    info.count[1] = emf->nx[1];

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
 
 Field solver
 
 *********************************************************************************************/

/**
 * @brief Advance magnetic field using Yee scheme
 * 
 * @param emf 	EM fields
 * @param dt 	Time step
 */
void yee_b( t_emf *emf, const float dt )
{
	// these must not be unsigned because we access negative cell indexes
	int i,j;
	float dt_dx, dt_dy;
	
    float3* const restrict B = emf -> B;
    const float3* const restrict E = emf -> E;

	dt_dx = dt / emf->dx[0];
	dt_dy = dt / emf->dx[1];
	
	// Canonical implementation
	const int nrow = emf->nrow;
	for (j=-1; j<=emf->nx[1]; j++) {
		for (i=-1; i<=emf->nx[0]; i++) {
			B[ i + j*nrow ].x += ( - dt_dy * ( E[i+(j+1)*nrow].z - E[i+j*nrow].z) );  
			B[ i + j*nrow ].y += (   dt_dx * ( E[(i+1)+j*nrow].z - E[i+j*nrow].z) );  
			B[ i + j*nrow ].z += ( - dt_dx * ( E[(i+1)+j*nrow].y - E[i+j*nrow].y) + 
									 dt_dy * ( E[i+(j+1)*nrow].x - E[i+j*nrow].x) );  
		}
	}
}

/**
 * @brief Advance Electric field using Yee scheme
 * 
 * @param emf 		EM fields
 * @param current 	Electric current density
 * @param dt 		Time step
 */
void yee_e( t_emf *emf, const t_current *current, const float dt )
{
	// these must not be unsigned because we access negative cell indexes
	int i,j;
	float dt_dx, dt_dy;
	
	dt_dx = dt / emf->dx[0];
	dt_dy = dt / emf->dx[1];

    float3* const restrict E = emf -> E;
    const float3* const restrict B = emf -> B;
    const float3* const restrict J = current -> J;
	
	// Canonical implementation
	const int nrow_e = emf->nrow;
	const int nrow_j = current->nrow;
	
	for (j=0; j<=emf->nx[1]+1; j++) {
		for (i=0; i<=emf->nx[0]+1; i++) {
			E[i+j*nrow_e].x += ( + dt_dy * ( B[i+j*nrow_e].z - B[i+(j-1)*nrow_e].z) ) 
			                     - dt * J[i+j*nrow_j].x;  
			
			E[i+j*nrow_e].y += ( - dt_dx * ( B[i+j*nrow_e].z - B[(i-1)+j*nrow_e].z) ) 
								 - dt * J[i+j*nrow_j].y;  

			E[i+j*nrow_e].z += ( + dt_dx * ( B[i+j*nrow_e].y - B[(i-1)+j*nrow_e].y) - 
								   dt_dy * ( B[i+j*nrow_e].x - B[i+(j-1)*nrow_e].x) ) 
									  - dt * J[i+j*nrow_j].z;  
 
		}
	}
}


/**
 * @brief Updates guard cell values
 *
 * When using periodic boundaries copies the lower cells to the upper guard
 * cells and vice-versa.
 * 
 * @param emf 
 */
void emf_update_gc( t_emf *emf )
{
	const int nrow = emf->nrow;

    float3* const restrict E = emf -> E;
    float3* const restrict B = emf -> B;

	// For moving window don't update x boundaries
	if ( ! emf -> moving_window ) {
		// x
		for (int j = -emf->gc[1][0]; j < emf->nx[1] + emf->gc[1][1]; j++) {
			
			// lower
			for (int i=-emf->gc[0][0]; i<0; i++) {
				E[ i + j*nrow ].x = E[ emf->nx[0] + i + j*nrow ].x;
				E[ i + j*nrow ].y = E[ emf->nx[0] + i + j*nrow ].y;
				E[ i + j*nrow ].z = E[ emf->nx[0] + i + j*nrow ].z;

				B[ i + j*nrow ].x = B[ emf->nx[0] + i + j*nrow ].x;
				B[ i + j*nrow ].y = B[ emf->nx[0] + i + j*nrow ].y;
				B[ i + j*nrow ].z = B[ emf->nx[0] + i + j*nrow ].z;
			}

			// upper
			for (int i=0; i<emf->gc[0][1]; i++) {
				E[ emf->nx[0] + i + j*nrow ].x = E[ i  + j*nrow ].x;
				E[ emf->nx[0] + i + j*nrow ].y = E[ i  + j*nrow ].y;
				E[ emf->nx[0] + i + j*nrow ].z = E[ i  + j*nrow ].z;
				
				B[ emf->nx[0] + i + j*nrow ].x = B[ i  + j*nrow ].x;
				B[ emf->nx[0] + i + j*nrow ].y = B[ i  + j*nrow ].y;
				B[ emf->nx[0] + i + j*nrow ].z = B[ i  + j*nrow ].z;
			}
			
		}
	}

	// y
	for (int i = -emf->gc[0][0]; i < emf->nx[0]+emf->gc[0][1]; i++) {
		
		// lower
		for (int j=-emf->gc[1][0]; j<0; j++) {
			E[ i + j*nrow ].x = E[ i + (emf->nx[1]+j)*nrow ].x;
			E[ i + j*nrow ].y = E[ i + (emf->nx[1]+j)*nrow ].y;
			E[ i + j*nrow ].z = E[ i + (emf->nx[1]+j)*nrow ].z;
			
			B[ i + j*nrow ].x = B[ i + (emf->nx[1]+j)*nrow ].x;
			B[ i + j*nrow ].y = B[ i + (emf->nx[1]+j)*nrow ].y;
			B[ i + j*nrow ].z = B[ i + (emf->nx[1]+j)*nrow ].z;
		}
		
		// upper
		for (int j=0; j<emf->gc[1][1]; j++) {
			E[ i + (emf->nx[1]+j)*nrow ].x = E[ i + j*nrow ].x;
			E[ i + (emf->nx[1]+j)*nrow ].y = E[ i + j*nrow ].y;
			E[ i + (emf->nx[1]+j)*nrow ].z = E[ i + j*nrow ].z;
			
			B[ i + (emf->nx[1]+j)*nrow ].x = B[ i + j*nrow ].x;
			B[ i + (emf->nx[1]+j)*nrow ].y = B[ i + j*nrow ].y;
			B[ i + (emf->nx[1]+j)*nrow ].z = B[ i + j*nrow ].z;
		}
		
	}
	
}

/**
 * @brief Move simulation window
 * 
 * When using a moving simulation window checks if a window move is due
 * at the current iteration and if so shifts left the data and zeroes
 * rightmost cells.
 * 
 * @param emf 
 */
void emf_move_window( t_emf *emf ){

	if ( ( emf -> iter * emf -> dt ) > emf->dx[0]*( emf -> n_move + 1 ) ) {
		const int nrow = emf->nrow;

	    float3* const restrict E = emf -> E;
	    float3* const restrict B = emf -> B;

	    const float3 zero_fld = {0.,0.,0.};

		// Shift data left 1 cell and zero rightmost cells
		for (int j = -emf->gc[1][0]; j < emf->nx[1] + emf->gc[1][1]; j++) {
			
			for (int i = -emf->gc[0][0]; i < emf->nx[0]+emf->gc[0][1] - 1; i++) {
				E[ i + j*nrow ] = E[ i + j*nrow + 1 ];
				B[ i + j*nrow ] = B[ i + j*nrow + 1 ];
			}

			for(int i = emf->nx[0]-1; i < emf->nx[0] + emf->gc[0][1]; i++) {
				E[ i + j*nrow ] = zero_fld;
				B[ i + j*nrow ] = zero_fld;				
			}
		}

		// Increase moving window counter
		emf -> n_move++;
	}
}

/**
 * @brief Advance EM fields 1 timestep
 * 
 * Fields are advanced in time using a FDTD algorith. The routine will also:
 * 1. Update guard cell values / apply boundary conditions
 * 2. Update "particle" fields if using external fields
 * 3. Move simulation window 
 * 
 * @param emf 		EM fields
 * @param current 	Electric current density
 */
void emf_advance( t_emf *emf, const t_current *current )
{
	uint64_t t0 = timer_ticks();
	const float dt = emf->dt;
	
	// Advance EM field using Yee algorithm modified for having E and B time centered
	yee_b( emf, dt/2.0f );
	
	yee_e( emf, current, dt );

	yee_b( emf, dt/2.0f );
	
	// Update guard cells with new values
	emf_update_gc( emf );

	// Update contribuition of external fields if necessary
	emf_update_part_fld( emf );

	// Advance internal iteration number
    emf -> iter += 1;

    // Move simulation window if needed
    if ( emf -> moving_window ) {
    	emf_move_window( emf );
    }
	
    // Update timing information
	_emf_time += timer_interval_seconds(t0, timer_ticks());
}

/**
 * @brief Calculate total EM field energy
 * 
 * Energy is calculated independently for each field component and is
 * returned as a 6 element vector for each of the E field components
 * [0..2] and B field components [3..5]. The energy is normalized to 
 * the cell size.
 * 
 * @param[in] emf EM fields
 * @param[out] energy Energy values vector
 */
void emf_get_energy( const t_emf *emf, double energy[] )
{
    float3* const restrict E = emf -> E;
    float3* const restrict B = emf -> B;
    const int nrow = emf -> nrow;

	for( int i = 0; i<6; i++) energy[i] = 0;

	for( int j = 0; j < emf -> nx[1]; j ++ ) {
		for( int i = 0; i < emf -> nx[0]; i ++ ) {
			energy[0] += E[i + j*nrow].x * E[i + j*nrow].x;
			energy[1] += E[i + j*nrow].y * E[i + j*nrow].y;
			energy[2] += E[i + j*nrow].z * E[i + j*nrow].z;
			energy[3] += B[i + j*nrow].x * B[i + j*nrow].x;
			energy[4] += B[i + j*nrow].y * B[i + j*nrow].y;
			energy[5] += B[i + j*nrow].z * B[i + j*nrow].z;
		}
	}

	for( int i = 0; i<6; i++) energy[i] *= 0.5 * emf -> dx[0] * emf -> dx[1];

}

/*********************************************************************************************

External Fields

 *********************************************************************************************/

/**
 * @brief Sets the external fields to be used for the simulation
 * 
 * @param emf 		EM field
 * @param ext_fld 	External fields
 */
void emf_set_ext_fld( t_emf* const emf, t_emf_ext_fld* ext_fld ) {

	emf -> ext_fld.E_type = ext_fld -> E_type;

	if ( emf -> ext_fld.E_type == EMF_FLD_TYPE_NONE ) {
		// Particle fields just point to the self-consistent fields
		emf -> E_part = emf -> E;
		emf -> ext_fld.E_part_buf = NULL;
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
        size_t size = (emf->gc[0][0] + emf->nx[0] + emf->gc[0][1]) * 
            (emf->gc[1][0] + emf->nx[1] + emf->gc[1][1]) * 
            sizeof( float3 );
	
		emf->ext_fld.E_part_buf = malloc( size );
        emf->E_part = emf->ext_fld.E_part_buf + emf->gc[0][0] + emf->gc[1][0] * emf->nrow;
	}

	emf -> ext_fld.B_type = ext_fld -> B_type;

	if ( emf -> ext_fld.B_type == EMF_FLD_TYPE_NONE ) {
		// Particle fields just point to the self-consistent fields
		emf -> B_part = emf -> B;
		emf -> ext_fld.B_part_buf = NULL;
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
        size_t size = (emf->gc[0][0] + emf->nx[0] + emf->gc[0][1]) * 
            (emf->gc[1][0] + emf->nx[1] + emf->gc[1][1]) * 
            sizeof( float3 );
	
		emf->ext_fld.B_part_buf = malloc( size );
        emf->B_part = emf->ext_fld.B_part_buf + emf->gc[0][0] + emf->gc[1][0] * emf->nrow;
	}

    // Initialize values on E/B_part grids
    emf_update_part_fld( emf );

}

/**
 * @brief Updates field values seen by particles with externally imposed fields
 * 
 * @param emf 	EM fields
 */
void emf_update_part_fld( t_emf* const emf ) {

    // Restrict pointers to E_part
    float3* const restrict E_part = emf->E_part;

    switch (emf->ext_fld.E_type)
    {
    case EMF_FLD_TYPE_UNIFORM: {
	    const int nrow_e = emf->nrow;
        for (int j=-emf->gc[1][0]; j<emf->nx[1]+emf->gc[1][1]; j++) {
            for (int i=-emf->gc[0][0]; i<emf->nx[0]+emf->gc[0][1]; i++) {
                float3 e = emf -> E[i+j*nrow_e];
                e.x += emf->ext_fld.E_0.x;
                e.y += emf->ext_fld.E_0.y;
                e.z += emf->ext_fld.E_0.z;
                E_part[i+j*nrow_e] = e;
            }
        }
        break; }
    case EMF_FLD_TYPE_CUSTOM: {
	    const int nrow_e = emf->nrow;
	
        for (int j=-emf->gc[1][0]; j<emf->nx[1]+emf->gc[1][1]; j++) {
            for (int i=-emf->gc[0][0]; i<emf->nx[0]+emf->gc[0][1]; i++) {
                float3 ext_E = (*emf->ext_fld.E_custom)
                    (i,emf->dx[0],j,emf->dx[1],emf->ext_fld.E_custom_data);

                float3 e = emf -> E[i+j*nrow_e];
                e.x += ext_E.x;
                e.y += ext_E.y;
                e.z += ext_E.z;
                E_part[i+j*nrow_e] = e;
            }
        }
        break; }
    case EMF_FLD_TYPE_NONE:
        break;
    }

    // Restrict pointers to B_part
    float3* const restrict B_part = emf->B_part;

    switch (emf->ext_fld.B_type)
    {
    case EMF_FLD_TYPE_UNIFORM: {
	    const int nrow_b = emf->nrow;
        for (int j=-emf->gc[1][0]; j<emf->nx[1]+emf->gc[1][1]; j++) {
            for (int i=-emf->gc[0][0]; i<emf->nx[0]+emf->gc[0][1]; i++) {
                float3 b = emf -> B[i+j*nrow_b];
                b.x += emf->ext_fld.B_0.x;
                b.y += emf->ext_fld.B_0.y;
                b.z += emf->ext_fld.B_0.z;
                B_part[i+j*nrow_b] = b;
            }
        }
        break; }
    case EMF_FLD_TYPE_CUSTOM: {
	    const int nrow_b = emf->nrow;
	
        for (int j=-emf->gc[1][0]; j<emf->nx[1]+emf->gc[1][1]; j++) {
            for (int i=-emf->gc[0][0]; i<emf->nx[0]+emf->gc[0][1]; i++) {
                float3 ext_B = (*emf->ext_fld.B_custom)
                    (i,emf->dx[0],j,emf->dx[1],emf->ext_fld.B_custom_data);

                float3 b = emf -> B[i+j*nrow_b];
                b.x += ext_B.x;
                b.y += ext_B.y;
                b.z += ext_B.z;
                B_part[i+j*nrow_b] = b;
            }
        }
        break; }
    case EMF_FLD_TYPE_NONE:
        break;
    }

}

/**
 * @brief Initialize EMF field values
 * 
 * @param emf       EM field object
 * @param init_fld  Initial field parameters
 */
void emf_init_fld( t_emf* const emf, t_emf_init_fld* init_fld )
{
    if ( emf -> iter != 0 ) {
        fprintf(stderr, "emf_init_fld should only be called at initialization, aborting...\n" );
        exit(-1);
    }

    float3* const restrict E = emf->E;
    float3* const restrict B = emf->B;
    const int stride = emf->nrow;

    switch ( init_fld -> E_type )
    {
    case EMF_FLD_TYPE_NONE:
        break;

    case EMF_FLD_TYPE_UNIFORM:
        for (int j=-emf->gc[1][0]; j<emf->nx[1]+emf->gc[1][1]; j++) {
            for (int i=-emf->gc[0][0]; i<emf->nx[0]+emf->gc[0][1]; i++) {
                E[ j*stride +  i ] = init_fld -> E_0;
            }
        }
        break;

    case EMF_FLD_TYPE_CUSTOM:
        for (int j=-emf->gc[1][0]; j<emf->nx[1]+emf->gc[1][1]; j++) {
            for (int i=-emf->gc[0][0]; i<emf->nx[0]+emf->gc[0][1]; i++) {
                float3 init_E = (init_fld->E_custom)
                    (i,emf->dx[0],j,emf->dx[1], init_fld->E_custom_data);
                E[ j*stride +  i ] = init_E;
            }
        }
        break;
    }    

    switch ( init_fld -> B_type )
    {
    case EMF_FLD_TYPE_NONE:
        break;

    case EMF_FLD_TYPE_UNIFORM:
        for (int j=-emf->gc[1][0]; j<emf->nx[1]+emf->gc[1][1]; j++) {
            for (int i=-emf->gc[0][0]; i<emf->nx[0]+emf->gc[0][1]; i++) {
                B[ j*stride +  i ] = init_fld -> B_0;
            }
        }
        break;

    case EMF_FLD_TYPE_CUSTOM:
        for (int j=-emf->gc[1][0]; j<emf->nx[1]+emf->gc[1][1]; j++) {
            for (int i=-emf->gc[0][0]; i<emf->nx[0]+emf->gc[0][1]; i++) {
                float3 init_B = (init_fld->B_custom)
                    (i,emf->dx[0],j,emf->dx[1], init_fld->B_custom_data);
                B[ j*stride +  i ] = init_B;
            }
        }
        break;
    }
}
