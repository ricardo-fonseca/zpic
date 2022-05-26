/**
 * @file current.c
 * @author Ricardo Fonseca
 * @brief Electric current density
 * @version 0.2
 * @date 2022-02-04
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include "current.h"

#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "zdf.h"

void current_smooth( t_current* const current );

/**
 * @brief Initializes Electric current density object
 * 
 * @param current   Electric current density
 * @param nx        Number of grid cells [x,y]
 * @param box       Physical box size [x,y]
 * @param dt        Simulation time step
 */
void current_new( t_current *current, int nx[], float box[], float dt )
{   
	// Number of guard cells for linear interpolation
	int gc[2][2] = {{1,2},
		            {1,2}}; 
	
	// Allocate global array
	size_t size;
	
	size = (gc[0][0] + nx[0] + gc[0][1]) * (gc[1][0] + nx[1] + gc[1][1]);
	
	current->J_buf = malloc( size * sizeof( float3 ) );
	assert( current->J_buf );

	// store nx and gc values
	for(int i = 0; i<2; i++){
		current->nx[i] = nx[i];
		current->gc[i][0] = gc[i][0];
		current->gc[i][1] = gc[i][1];
	}
	current -> nrow = gc[0][0] + nx[0] + gc[0][1];
	
	// Make J point to cell [0][0]
	current->J = current->J_buf + gc[0][0] + gc[1][0] * current->nrow;
	
	// Set cell sizes and box limits
	for(int i = 0; i<2; i++){
		current -> box[i] = box[i];
		current -> dx[i]  = box[i] / nx[i];
	}

	// Clear smoothing options
	current -> smooth = (t_smooth) {
		.xtype = NONE,
		.ytype = NONE,
		.xlevel = 0,
		.ylevel = 0
	};

	// Initialize time information
	current -> iter = 0;
	current -> dt = dt;

	current -> moving_window = 0;

    // Zero initial current
    // This is only relevant for diagnostics, current is always zeroed before deposition
	current_zero( current );
	
}

/**
 * @brief Frees dynamic memory from electric current density
 * 
 * @param current   Electric current density
 */
void current_delete( t_current *current )
{
	free( current->J_buf );
	
	current->J_buf = NULL;
}

/**
 * @brief Sets all electric current density values to zero
 * 
 * @param current   Electric current density
 */
void current_zero( t_current *current )
{
	// zero fields
	size_t size;
	
	size = (current->gc[0][0] + current->nx[0] + current->gc[0][1]) * 
	       (current->gc[1][0] + current->nx[1] + current->gc[1][1]) * sizeof( float3 );
	memset( current->J_buf, 0, size );
	
}

/**
 * @brief Updates guard cell values
 * 
 * When using periodic boundaries the electric current that was added to
 * the upper guard cells will be added to the corresponding lower grid
 * cells, and the values then copied to the upper grid cells
 * 
 * @param current Electric current density
 */
void current_update_gc( t_current *current )
{
	const int nrow = current->nrow;
	float3* restrict const J = current -> J;
	
	// x
	if ( ! current -> moving_window ) {
		for (int j = -current->gc[1][0]; j < current->nx[1] + current->gc[1][1]; j++) {
			
			// lower - add the values from upper boundary ( both gc and inside box )
			for (int i=-current->gc[0][0]; i<current->gc[0][1]; i++) {
				J[ i + j*nrow ].x += J[ current->nx[0] + i + j*nrow ].x;
				J[ i + j*nrow ].y += J[ current->nx[0] + i + j*nrow ].y;
				J[ i + j*nrow ].z += J[ current->nx[0] + i + j*nrow ].z;
			}
			
			// upper - just copy the values from the lower boundary 
			for (int i=-current->gc[0][0]; i<current->gc[0][1]; i++) {
				J[ current->nx[0] + i + j*nrow ] = J[ i  + j*nrow ];
			}
			
		}
	}
	
	// y
	for (int i = -current->gc[0][0]; i < current->nx[0]+current->gc[0][1]; i++) {
		
		// lower - add the values from upper boundary ( both gc and inside box )
		for (int j=-current->gc[1][0]; j<current->gc[1][1]; j++) {
			J[ i + j*nrow ].x += J[ i + (current->nx[1]+j)*nrow ].x;
			J[ i + j*nrow ].y += J[ i + (current->nx[1]+j)*nrow ].y;
			J[ i + j*nrow ].z += J[ i + (current->nx[1]+j)*nrow ].z;
		}
		
		// upper - just copy the values from the lower boundary 
		for (int j=-current->gc[1][0]; j<current->gc[1][1]; j++) {
			J[ i + (current->nx[1]+j)*nrow ] = J[ i + j*nrow ];
		}
		
	}
	
}

/**
 * @brief Advances electric current density 1 time step
 * 
 * The routine will:
 * 1. Update the guard cells
 * 2. Apply digitial filtering (if configured)
 * 3. Advance iteration number
 * 
 * @param current Electric current density
 */
void current_update( t_current *current )
{
    
    // Boundary conditions / guard cells
    current_update_gc( current );

    // Smoothing
    current_smooth( current );

    // Advance iteration number
    current -> iter++;
    
}

/**
 * @brief Saves electric current density diagnostic information to disk
 * 
 * Saves the selected current density component to disk in directory
 * "CURRENT". Guard cell values are discarded.
 * 
 * @param current Electric current object
 * @param jc Current component to save, must be one of {0,1,2}
 */
void current_report( const t_current *current, const int jc )
{
	if ( jc < 0 || jc > 2 ) {
		fprintf(stderr, "(*error*) Invalid current component (jc) selected, returning\n");
		return;
	}
		
	// Pack the information
	float * buf = malloc( current->nx[0]*current->nx[1]*sizeof(float) );
    float * p = buf;
	float3 * f = current->J;
	
	switch (jc) {
		case 0:
			for( int j = 0; j < current->nx[1]; j++) {
				for ( int i = 0; i < current->nx[0]; i++ ) {
					p[i] = f[i].x;
				}
				p += current->nx[0];
				f += current->nrow;
			}
			break;
		case 1:
			for( int j = 0; j < current->nx[1]; j++) {
				for ( int i = 0; i < current->nx[0]; i++ ) {
					p[i] = f[i].y;
				}
				p += current->nx[0];
				f += current->nrow;
			}
			break;
		case 2:
			for( int j = 0; j < current->nx[1]; j++) {
				for ( int i = 0; i < current->nx[0]; i++ ) {
					p[i] = f[i].z;
				}
				p += current->nx[0];
				f += current->nrow;
			}
			break;
	}

	char vfname[16];	// Dataset name
	char vflabel[16];	// Dataset label (for plots)

    snprintf( vfname, 3, "J%1u", jc );
	char comp[] = {'x','y','z'};
    snprintf(vflabel,4,"J_%c",comp[jc]);
	
    t_zdf_grid_axis axis[2];
    axis[0] = (t_zdf_grid_axis) {
    	.min = 0.0,
    	.max = current->box[0],
		.name = "x",
    	.label = "x",
    	.units = "c/\\omega_p"
    };

    axis[1] = (t_zdf_grid_axis) {
    	.min = 0.0,
    	.max = current->box[1],
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

    info.count[0] = current->nx[0];
    info.count[1] = current->nx[1];

    t_zdf_iteration iter = {
		.name = "ITERATION",
    	.n = current->iter,
    	.t = current -> iter * current -> dt,
    	.time_units = "1/\\omega_p"
    };

	zdf_save_grid( (void *) buf, zdf_float32, &info, &iter, "CURRENT" );
	
	// free local data
	free( buf );
	
}

/**
 * @brief Gets the value of the compensator kernel for an n pass binomial kernel
 * 
 * This kernel eliminates the $k^2$ dependency of the transfer function
 * near $k = 0$. The resulting kernel will be in the form [a,b,a], with
 * the values of a and b being determined by this function. The result
 * is normalized.
 * 
 * @param n Number of binomial passes
 * @param sa a value of the compensator kernel
 * @param sb b value of the compensator kernel
 */
void get_smooth_comp( int n, float* sa, float* sb) {
	float a = -1;
	float b = (4.0 + 2.0*n)/n;
	float total = 2*a + b;

	*sa = a / total;
	*sb = b / total;
}

/**
 * @brief Applies a 3 point kernel convolution along x
 * 
 * The kernel has the form [a,b,a]. The routine accounts for periodic
 * boundaries.
 * 
 * @param current 
 * @param sa kernel a value
 * @param sb kernel b value
 */
void kernel_x( t_current* const current, const float sa, const float sb )
{
	float3* restrict const J = current -> J;
	const int nrow = current->nrow;

	for( int j = 0; j < current -> nx[1]; j++) {
		
		int idx = j * nrow;

		float3 fl = J[idx-1];
		float3 f0 = J[idx  ];

		for( int i = 0; i < current -> nx[0]; i++) {

			float3 fu = J[idx + i + 1];

			float3 fs;

			fs.x = sa * fl.x + sb * f0.x + sa * fu.x;
			fs.y = sa * fl.y + sb * f0.y + sa * fu.y;
			fs.z = sa * fl.z + sb * f0.z + sa * fu.z;

			J[idx + i] = fs;

			fl = f0;
			f0 = fu;

		}

		// Update x boundaries unless we are using a moving window
		if ( ! current -> moving_window ) {
			for(int i = -current->gc[0][0]; i<0; i++) 
				J[ idx + i ] = J[ idx + current->nx[0] + i ];

			for (int i=0; i<current->gc[0][1]; i++)
				J[ idx + current->nx[0] + i ] = J[ idx + i ];
		}
	}
}

/**
 * @brief Applies a 3 point kernel convolution along y
 * 
 * The kernel has the form [a,b,a]. The routine accounts for periodic
 * boundaries.
 * 
 * @param current 
 * @param sa kernel a value
 * @param sb kernel b value
 */
void kernel_y( t_current* const current, const float sa, const float sb ){

	float3 flbuf[ current -> nx[0] ];
	float3* restrict const J = current -> J;
	const int nrow = current->nrow;

	// buffer lower row
	for( int i = 0; i < current -> nx[0]; i++) {
		flbuf[i] = J[i - nrow];
	}


	for( int j = 0; j < current -> nx[1]; j++) {
		
		int idx = j * nrow;

		for( int i = 0; i < current -> nx[0]; i++) {

			// Get lower, central and upper values
			float3 fl = flbuf[i];
			float3 f0 = J[ idx + i ];
			float3 fu = J[ idx + i + nrow ];

			// Store the value that will be overritten for use in the next row
			flbuf[i] = f0;

			// Convolution with kernel
			float3 fs;
			fs.x = sa * fl.x + sb * f0.x + sa * fu.x;
			fs.y = sa * fl.y + sb * f0.y + sa * fu.y;
			fs.z = sa * fl.z + sb * f0.z + sa * fu.z;

			// Store result
			J[idx + i] = fs;
		}
	}

	// Update y boundaries

	// Grid is always periodic along y
	for (int i = -current->gc[0][0]; i < current->nx[0]+current->gc[0][1]; i++) {
		for (int j=-current->gc[1][0]; j<0; j++) 
			J[ i + j*nrow ] = J[ i + (current->nx[1]+j)*nrow ];
		for (int j=0; j<current->gc[1][1]; j++)
			J[ i + (current->nx[1]+j)*nrow ] = J[ i + j*nrow ];
	}

}

/**
 * @brief Applies digital filtering to the current density
 * 
 * Filtering is applied through a sequence of 3 point kernel convolutions.
 * The routine will apply a binomial kernel ([1,2,1]) n times, followed by
 * an optional compensator kernel. X and Y directions are treated
 * separately
 * 
 * Filtering parameters are set by the `current -> smooth` variable.
 * 
 * @param current Electric current density
 */
void current_smooth( t_current* const current )
{
	// filter kernel [sa, sb, sa]
	float sa, sb;

    // x-direction filtering
    if ( current -> smooth.xtype != NONE ) {
    	// binomial filter
    	for( int i = 0; i < current -> smooth.xlevel; i++) {
    		kernel_x( current, 0.25, 0.5 );
    	}

    	// Compensator
    	if ( current -> smooth.xtype == COMPENSATED ) {
    		get_smooth_comp( current -> smooth.xlevel, &sa, &sb );
    		kernel_x( current, sa, sb );
    	}
    }

    // y-direction filtering
    if ( current -> smooth.ytype != NONE ) {
    	// binomial filter
    	for( int i = 0; i < current -> smooth.xlevel; i++) {
    		kernel_y( current, 0.25, 0.5 );
    	}

    	// Compensator
    	if ( current -> smooth.ytype == COMPENSATED ) {
    		get_smooth_comp( current -> smooth.ylevel, &sa, &sb );
    		kernel_y( current, sa, sb );
    	}
    }
}

