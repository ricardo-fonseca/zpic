/*
 *  particles.c
 *  zpic
 *
 *  Created by Ricardo Fonseca on 11/8/10.
 *  Copyright 2010 Centro de Física dos Plasmas. All rights reserved.
 *
 */

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

#include "grid2d.h"

#include "particles.h"

#include "random.h"
#include "emf.h"
#include "current.h"

#include "zdf.h"
#include "timer.h"

static double _spec_time = 0.0;
static double _spec_npush = 0.0;

void spec_sort( t_species *spec );

/**
 * Returns the total time spent pushing particles (includes boundaries and moving window)
 * @return  Total time in seconds
 */
double spec_time( void )
{
	return _spec_time;
}

/**
 * Returns the performance achieved by the code (push time)
 * @return  Performance in seconds per particle
 */
double spec_perf( void )
{
	return (_spec_npush > 0 )? _spec_time / _spec_npush: 0.0;
}

/*********************************************************************************************
 
 Initialization
 
 *********************************************************************************************/

void spec_set_u( t_species* spec, const int start, const int end )
{
	int i;    

	for (i = start; i <= end; i++) {
		spec->part[i].ux = spec -> ufl[0] + spec -> uth[0] * rand_norm(); 
		spec->part[i].uy = spec -> ufl[1] + spec -> uth[1] * rand_norm(); 
		spec->part[i].uz = spec -> ufl[2] + spec -> uth[2] * rand_norm(); 
	}

}	

void spec_set_x( t_species* spec, const int range[][2] )
{

	int i, j, k, ip;
	
	float* poscell;
	float start, end;
	
	// Calculate particle positions inside the cell
	const int npc = spec->ppc[0]*spec->ppc[1];
	
	poscell = malloc( 2 * npc * sizeof( float ) );
	ip = 0;
	for (j =0; j<spec->ppc[1]; j++) {
		for (i=0; i<spec->ppc[0]; i++) {
			poscell[ip]   = (1 + 2*i - spec->ppc[0]) / (2.0*spec->ppc[0]);
			poscell[ip+1] = (1 + 2*j - spec->ppc[1]) / (2.0*spec->ppc[1]);
			ip+=2;
		}
	}

	ip = spec -> np;
	
	// Set position of particles in the specified grid range according to the density profile
	switch ( spec -> density.type ) {
	case STEP: // Step like density profile
		
		// Get edge position normalized to cell size;
		start = spec -> density.start / spec -> dx[0];

		for (j = range[1][0]; j <= range[1][1]; j++) {
			for (i = range[0][0]; i <= range[0][1]; i++) {

				for (k=0; k<npc; k++) {
					if ( i + poscell[2*k] > start ) {
						spec->part[ip].ix = i;
						spec->part[ip].iy = j;
						spec->part[ip].x = poscell[2*k];
						spec->part[ip].y = poscell[2*k+1];
						ip++;
					}
				}
			}
		}
		break;

	case SLAB: // Slab like density profile
		
		// Get edges positions normalized to cell size;
		start = spec -> density.start / spec -> dx[0] - 0.5;
		end   = spec -> density.end / spec -> dx[0] - 0.5;

		for (j = range[1][0]; j <= range[1][1]; j++) {
			for (i = range[0][0]; i <= range[0][1]; i++) {

				for (k=0; k<npc; k++) {
					if ( i + poscell[2*k] > start &&  i + poscell[2*k] < end ) {
						spec->part[ip].ix = i;
						spec->part[ip].iy = j;
						spec->part[ip].x = poscell[2*k];
						spec->part[ip].y = poscell[2*k+1];
						ip++;
					}
				}
			}
		}
		break;

	default: // Uniform density
		for (j = range[1][0]; j <= range[1][1]; j++) {
			for (i = range[0][0]; i <= range[0][1]; i++) {

				for (k=0; k<npc; k++) {
					spec->part[ip].ix = i;
					spec->part[ip].iy = j;
					spec->part[ip].x = poscell[2*k];
					spec->part[ip].y = poscell[2*k+1];
					ip++;
				}
			}
		}
	}
	
	spec -> np = ip;
	
	free(poscell);
	
}

void spec_inject_particles( t_species* spec, const int range[][2] )
{
	int start = spec -> np;

	// Get maximum number of particles to inject
	int np_inj = ( range[0][1] - range[0][0] + 1 ) * ( range[1][1] - range[1][0] + 1 ) * 
	             spec -> ppc[0] * spec -> ppc[1];

	// Check if buffer is large enough and if not reallocate
	if ( spec -> np + np_inj > spec -> np_max ) {
        spec -> np_max = (( spec -> np_max + np_inj )/1024 + 1) * 1024;
		spec -> part = realloc( (void*) spec -> part, spec -> np_max * sizeof(t_part) );
	}

	// Set particle positions
	spec_set_x( spec, range );
	
	// Set momentum of injected particles
	spec_set_u( spec, start, spec -> np - 1 );

}

void spec_new( t_species* spec, char name[], const float m_q, const int ppc[], 
			  const float *ufl, const float * uth,
			  const int nx[], float box[], const float dt, t_density* density )
{

	int i, npc;
	
	// Species name
	strncpy( spec -> name, name, MAX_SPNAME_LEN );
	
	npc = 1;
	// Store species data
	for (i=0; i<2; i++) {
		spec->nx[i] = nx[i];
		spec->ppc[i] = ppc[i];
		npc *= ppc[i];
		
		spec->box[i] = box[i];
		spec->dx[i] = box[i] / nx[i];
	}
	
	spec -> m_q = m_q;
	spec -> q = copysign( 1.0f, m_q ) / npc;

	spec -> dt = dt;
	
	// Initialize particle buffer
	spec->np_max = 0;
	spec->part = NULL;
	
	
	// Initialize density profile
	if ( density ) {
		spec -> density = *density;
		if ( spec -> density.n == 0. ) spec -> density.n = 1.0;
	} else {
		// Default values
		spec -> density = (t_density) { .type = UNIFORM, .n = 1.0 };
	}

	// Density multiplier
	spec ->q *= fabsf( spec -> density.n );

	// Initialize temperature profile
	if ( ufl ) {
		for(i=0; i<3; i++) spec -> ufl[i] = ufl[i];
	} else {
		for(i=0; i<3; i++) spec -> ufl[i] = 0;
	}

	if ( uth ) {
		for(i=0; i<3; i++) spec -> uth[i] = uth[i];
	} else {
		for(i=0; i<3; i++) spec -> uth[i] = 0;
	}

	// Reset iteration number
	spec -> iter = 0;

    // Inject initial particle distribution
    spec -> np = 0;
    
    const int range[][2] = {{0, nx[0]-1},
                            {0, nx[1]-1}};

    spec_inject_particles( spec, range );

}

void spec_delete( t_species* spec )
{
	free(spec->part);
	spec->np = -1;
}


/*********************************************************************************************
 
 Charge / Current deposition
 
 *********************************************************************************************/

int ltrim( float x )
{
	return (( x >= 0.5f )?1:0) - (( x < -0.5f )?1:0);
}


void deposit_charge( t_scalar_grid2d * rho, const t_part* restrict const part, const float q )
{
	float s0x, s1x;
	float s0y, s1y;

	const int nrow = rho -> nrow;
	const int idx  = part->ix + part->iy * nrow;
	
	s0x = 0.5f - part->x;
	s1x = 0.5f + part->x; 

	s0y = 0.5f - part->y;
	s1y = 0.5f + part->y; 

	rho->s[ idx            ] += s0y * s0x * q;
	rho->s[ idx        + 1 ] += s0y * s1x * q;
	rho->s[ idx + nrow     ] += s1y * s0x * q;
	rho->s[ idx + nrow + 1 ] += s1y * s1x * q;

}

void deposit_current( t_vfld_grid2d* J, const t_part* restrict const part, const float q, const float rg, 
	const float dx, const float dy )
{
	int i, di;
	int j, dj;
	float x, y;
	float s0x, s1x, s0y, s1y; 
	float jx, jy, jz;

	const int nrow = J -> nrow;
	
	// Find position time centered with velocity
	i = part->ix;
	x = part->x + 0.5f * dx;
	di = ltrim(x);
	i += di;
	x -= di;

	j = part->iy;
	y = part->y + 0.5f * dy;
	dj = ltrim(y);
	j += dj;
	y -= dj;
	
	s0x = 0.5f - x;
	s1x = 0.5f + x; 

	s0y = 0.5f - y;
	s1y = 0.5f + y; 

	jx = q * part -> ux * rg;
	jy = q * part -> uy * rg;
	jz = q * part -> uz * rg;

	int idx = i + j * nrow;

	J->x[idx           ] += s0y * s0x * jx;
	J->y[idx           ] += s0y * s0x * jy;
	J->z[idx           ] += s0y * s0x * jz;

	J->x[idx        + 1] += s0y * s1x * jx;
	J->y[idx        + 1] += s0y * s1x * jy;
	J->z[idx        + 1] += s0y * s1x * jz;

	J->x[idx + nrow    ] += s1y * s0x * jx;
	J->y[idx + nrow    ] += s1y * s0x * jy;
	J->z[idx + nrow    ] += s1y * s0x * jz;

	J->x[idx + nrow + 1] += s1y * s1x * jx;
	J->y[idx + nrow + 1] += s1y * s1x * jy;
	J->z[idx + nrow + 1] += s1y * s1x * jz;

}

/*********************************************************************************************
 
 Sorting
 
 *********************************************************************************************/

void spec_sort( t_species* spec )
{
	int *idx, *npic;

	int ncell = spec->nx[0]*spec->nx[1];
	
	// Allocate index memory
	idx  = malloc(spec->np*sizeof(int));

	// Allocate temp. array with number of particles in cell
	npic = malloc( ncell * sizeof(int));
	memset( npic, 0, ncell * sizeof(int));

	// Generate sorted index
    int i;
	for (i=0; i<spec->np; i++) {
		idx[i] = spec->part[i].ix + spec->part[i].iy * spec->nx[0];
		npic[idx[i]]++;
	}
	
	int isum = 0, j;
	for (i=0; i<ncell; i++) {
		j = npic[i];
		npic[i] = isum;
		isum += j;
	}
	
	for (i=0; i< spec->np; i++) {
		j = idx[i];
		idx[i] = npic[j]++;
	}
	
	// free temp. array
	free(npic);

	// low mem
	for (i=0; i < spec->np; i++) {
		register t_part tmp;
		register int k;
		
		k = idx[i];
		while ( k > i ) {
			register int t;
			
			tmp = spec->part[k];
			spec->part[k] = spec->part[i];
			spec->part[i] = tmp;
			
			t = idx[k];
			idx[k] = -1;
			k = t;
		}
	}

	free(idx);

}


/*********************************************************************************************
 
 Particle advance
 
 *********************************************************************************************/


void interpolate_fld( t_vfld_grid2d * E, t_vfld_grid2d * B, 
	          const t_part* restrict const part, t_vfld* restrict const Ep, t_vfld* restrict const Bp )
{
	float s0x, s0y, s1x, s1y;

	const int nrow = E -> nrow;
	
	const int idx = part->ix + part->iy * nrow;
	
	s0x = 0.5f - part->x;
	s1x = 0.5f + part->x;

	s0y = 0.5f - part->y;
	s1y = 0.5f + part->y;

	Ep -> x = s0y * s0x * E->x[idx           ] + 
	          s0y * s1x * E->x[idx        + 1] + 
	          s1y * s0x * E->x[idx + nrow    ] + 
	          s1y * s1x * E->x[idx + nrow + 1];

	Ep -> y = s0y * s0x * E->y[idx           ] +
	          s0y * s1x * E->y[idx        + 1] + 
	          s1y * s0x * E->y[idx + nrow    ] + 
	          s1y * s1x * E->y[idx + nrow + 1];

	Ep -> z = s0y * s0x * E->z[idx           ] + 
	          s0y * s1x * E->z[idx        + 1] + 
	          s1y * s0x * E->z[idx + nrow    ] +
	          s1y * s1x * E->z[idx + nrow + 1];

	Bp -> x = s0y * s0x * B->x[idx           ] + 
	          s0y * s1x * B->x[idx        + 1] + 
	          s1y * s0x * B->x[idx + nrow    ] + 
	          s1y * s1x * B->x[idx + nrow + 1];

	Bp -> y = s0y * s0x * B->y[idx           ] + 
	          s0y * s1x * B->y[idx        + 1] + 
	          s1y * s0x * B->y[idx + nrow    ] + 
	          s1y * s1x * B->y[idx + nrow + 1];

	Bp -> z = s0y * s0x * B->z[idx           ] + 
	          s0y * s1x * B->z[idx        + 1] + 
	          s1y * s0x * B->z[idx + nrow    ] +
	          s1y * s1x * B->z[idx + nrow + 1];

}	


void spec_advance( t_species* spec, t_emf* emf, t_charge* charge, t_current* current )
{
	int i;
	
	uint64_t t0;
	t0 = timer_ticks();
	
	const float tem   = 0.5 * spec->dt/spec -> m_q;
	const float dt_dx = spec->dt / spec->dx[0]; 
	const float dt_dy = spec->dt / spec->dx[1]; 

	const int nx0 = spec -> nx[0];
	const int nx1 = spec -> nx[1];

	double energy = 0;

	// Advance particles
	for (i=0; i<spec->np; i++) {
				
		t_vfld Ep, Bp;
		float utx, uty, utz;
		float ux, uy, uz, u2;
		float gamma, rg, gtem, otsq;
		
		float x1, y1;
		
		int di, dj;
		float dx, dy; 

		// Load particle momenta
		ux = spec -> part[i].ux;
		uy = spec -> part[i].uy;
		uz = spec -> part[i].uz;

		// interpolate fields
		interpolate_fld( &emf -> E, &emf -> B, &spec -> part[i], &Ep, &Bp );
		
		// advance u using Boris scheme
		Ep.x *= tem;
		Ep.y *= tem;
		Ep.z *= tem;
		
		utx = ux + Ep.x;
		uty = uy + Ep.y;
		utz = uz + Ep.z;

		// Perform first half of the rotation

		// Get time centered gamma
		u2 = utx*utx + uty*uty + utz*utz;
		gamma = sqrtf( 1 + u2 );

		// Accumulate time centered energy
		energy += u2 / ( 1 + gamma );

		gtem = tem / gamma;
		
		Bp.x *= gtem;
		Bp.y *= gtem;
		Bp.z *= gtem;
		otsq = 2.0f / ( 1.0f + Bp.x*Bp.x + Bp.y*Bp.y + Bp.z*Bp.z );

		ux = utx + uty*Bp.z - utz*Bp.y;
		uy = uty + utz*Bp.x - utx*Bp.z;
		uz = utz + utx*Bp.y - uty*Bp.x;
		
		// Perform second half of the rotation
		
		Bp.x *= otsq;
		Bp.y *= otsq;
		Bp.z *= otsq;
		
		utx += uy*Bp.z - uz*Bp.y;
		uty += uz*Bp.x - ux*Bp.z;
		utz += ux*Bp.y - uy*Bp.x;
		
		// Perform second half of electric field acceleration
		ux = utx + Ep.x;
		uy = uty + Ep.y;
		uz = utz + Ep.z;
		
		// Store new momenta
		spec -> part[i].ux = ux;
		spec -> part[i].uy = uy;
		spec -> part[i].uz = uz;
		
		// push particle
		rg = 1.0f / sqrtf(1.0f + ux*ux + uy*uy + uz*uz);
				
		dx = dt_dx * rg * ux;
		dy = dt_dy * rg * uy;

		// Deposit current for the particle (t+1/2)
		deposit_current( &current -> J, &spec -> part[i], spec->q, rg, dx, dy );
		
		x1 = spec -> part[i].x + dx; 
		y1 = spec -> part[i].y + dy;
		
		di = ltrim(x1);
		dj = ltrim(y1);

		x1 -= di;
		y1 -= dj;
		
		
		// Store results
		spec -> part[i].x = x1;
		spec -> part[i].y = y1;
		spec -> part[i].ix += di;
		spec -> part[i].iy += dj;

		// Deposit charge for the particle (t+1)
		deposit_charge( &charge -> rho, &spec -> part[i], spec->q );

		
	}

	// Store energy
	spec -> energy = spec-> q * spec -> m_q * energy * spec -> dx[0] * spec -> dx[1];

	// Advance internal iteration number
    spec -> iter += 1;
    _spec_npush += spec -> np;

	// Use periodic boundaries in both directions
	for (i=0; i<spec->np; i++) {
		spec -> part[i].ix += (( spec -> part[i].ix < 0 ) ? nx0 : 0 ) - (( spec -> part[i].ix >= nx0 ) ? nx0 : 0);
		spec -> part[i].iy += (( spec -> part[i].iy < 0 ) ? nx1 : 0 ) - (( spec -> part[i].iy >= nx1 ) ? nx1 : 0);
	}
	
	// Sort species at every 16 time steps
	if ( ! (spec -> iter % 16) ) spec_sort( spec );
    	
	
	_spec_time += timer_interval_seconds( t0, timer_ticks() );
}


/*********************************************************************************************
 
 Charge Deposition
 
 *********************************************************************************************/


void spec_deposit_charge( const t_species* spec, float* charge )
{
	int i,j;
	
	// Charge array is expected to have 1 guard cell at the upper boundary
	const int nrow = spec -> nx[0] + 1;
	const float q = spec -> q;
	
	for (i=0; i<spec->np; i++) {
		int idx = spec->part[i].ix + nrow*spec->part[i].iy;

		float s0x = 0.5f - spec->part[i].x;
		float s1x = 0.5f + spec->part[i].x; 

		float s0y = 0.5f - spec->part[i].y;
		float s1y = 0.5f + spec->part[i].y; 
		
		charge[ idx            ] += s0y * s0x * q;
		charge[ idx + 1        ] += s0y * s1x * q;
		charge[ idx     + nrow ] += s1y * s0x * q;
		charge[ idx + 1 + nrow ] += s1y * s1x * q;
	}

	// Correct boundary values

	// x - Periodic boundaries
	for (j = 0; j < spec -> nx[1] + 1; j++) {
		charge[ 0 + j*nrow ] += charge[ spec -> nx[0] + j*nrow ];
	}
	
	// y - Periodic boundaries
	for (i = 0; i < spec->nx[0]+1; i++) {
		charge[ i + 0 ] += charge[ i + spec -> nx[1] * nrow ];
	}

}

/*********************************************************************************************
 
 Diagnostics
 
 *********************************************************************************************/

void spec_rep_particles( const t_species *spec )
{
	
	t_zdf_file part_file;

	int i;
	
	const char * quants[] = {
	    "x1","x2",
	    "u1","u2","u3"
	};

	const char * units[] = {
	    "c/\\omega_p", "c/\\omega_p",
	    "c","c","c"
	};

    t_zdf_iteration iter = {
    	.n = spec->iter,
    	.t = spec -> iter * spec -> dt,
    	.time_units = "1/\\omega_p"
    };

	// Allocate buffer for positions
	
	t_zdf_part_info info = {
		.name = (char *) spec -> name,
		.nquants = 5,
		.quants = (char **) quants,
		.units = (char **) units,
		.np = spec ->np
	};

	// Create file and add description
	zdf_part_file_open( &part_file, &info, &iter, "PARTICLES" );

	// Add positions and generalized velocities
	size_t size = ( spec -> np ) * sizeof( float );
	float* data = malloc( size );

	// x1
	for( i = 0; i < spec ->np; i++ )
		data[i] = ( spec -> part[i].ix + (spec -> part[i].x +0.5f ) ) * spec -> dx[0];
	zdf_part_file_add_quant( &part_file, quants[0], data, spec ->np );

	// x2
	for( i = 0; i < spec ->np; i++ )
		data[i] = ( spec -> part[i].iy + (spec -> part[i].y + 0.5f ) ) * spec -> dx[1];
	zdf_part_file_add_quant( &part_file, quants[1], data, spec ->np );

	// ux
	for( i = 0; i < spec ->np; i++ ) data[i] = spec -> part[i].ux;
	zdf_part_file_add_quant( &part_file, quants[2], data, spec ->np );

	// uy
	for( i = 0; i < spec ->np; i++ ) data[i] = spec -> part[i].uy;
	zdf_part_file_add_quant( &part_file, quants[3], data, spec ->np );

	// uz
	for( i = 0; i < spec ->np; i++ ) data[i] = spec -> part[i].uz;
	zdf_part_file_add_quant( &part_file, quants[4], data, spec ->np );

	free( data );

	zdf_close_file( &part_file );
}	


void spec_rep_charge( const t_species *spec )
{
	float *buf, *charge, *b, *c;
	size_t size;
	int i, j;
	
	// Add 1 guard cell to the upper boundary
	size = ( spec -> nx[0] + 1 ) * ( spec -> nx[1] + 1 ) * sizeof( float );
	charge = malloc( size );
	memset( charge, 0, size );
	
	// Deposit the charge
	spec_deposit_charge( spec, charge );
	
	// Compact the data to save the file (throw away guard cells)
	size = ( spec -> nx[0] ) * ( spec -> nx[1] );
	buf = malloc( size * sizeof( float ) );
	
	b = buf;
	c = charge;
	for( j = 0; j < spec->nx[1]; j++) {
		for ( i = 0; i < spec->nx[0]; i++ ) {
			b[i] = c[i];
		}
		b += spec->nx[0];
		c += spec->nx[0] + 1;
	}
	
	free( charge );

    t_zdf_grid_axis axis[2];
    axis[0] = (t_zdf_grid_axis) {
    	.min = 0.0,
    	.max = spec->box[0],
    	.label = "x_1",
    	.units = "c/\\omega_p"
    };

    axis[1] = (t_zdf_grid_axis) {
    	.min = 0.0,
    	.max = spec->box[1],
    	.label = "x_2",
    	.units = "c/\\omega_p"
    };

    t_zdf_grid_info info = {
    	.ndims = 2,
    	.label = "charge",
    	.units = "n_e",
    	.axis  = axis
    };

    info.nx[0] = spec->nx[0];
    info.nx[1] = spec->nx[1];

    t_zdf_iteration iter = {
    	.n = spec->iter,
    	.t = spec -> iter * spec -> dt,
    	.time_units = "1/\\omega_p"
    };

	zdf_save_grid( buf, &info, &iter, spec->name );	


	free( buf );
}	


void spec_pha_axis( const t_species *spec, int i0, int np, int quant, float *axis )
{
	int i;
	
	switch (quant) {
		case X1:
			for (i = 0; i < np; i++) 
				axis[i] = ( (spec -> part[i0+i].x + 0.5f) + spec -> part[i0+i].ix ) * spec -> dx[0];
			break;
		case X2:
			for (i = 0; i < np; i++) 
				axis[i] = ( (spec -> part[i0+i].y + 0.5f) + spec -> part[i0+i].iy ) * spec -> dx[1];
			break;
		case U1:
			for (i = 0; i < np; i++) 
				axis[i] = spec -> part[i0+i].ux;
			break;
		case U2:
			for (i = 0; i < np; i++) 
				axis[i] = spec -> part[i0+i].uy;
			break;
		case U3:
			for (i = 0; i < np; i++) 
				axis[i] = spec -> part[i0+i].uz;
			break;
	}
}

const char * spec_pha_axis_units( int quant ) {
	switch (quant) {
		case X1:
			return("c/\\omega_p");
			break;
		case U1:
		case U2:
		case U3:
			return("m_e c");
	}
	return("");
}

void spec_deposit_pha( const t_species *spec, const int rep_type,
			  const int pha_nx[], const float pha_range[][2], float* restrict buf )
{
	const int BUF_SIZE = 1024;
	float pha_x1[BUF_SIZE], pha_x2[BUF_SIZE];


	const int nrow = pha_nx[0];

	const int quant1 = rep_type & 0x000F;
	const int quant2 = (rep_type & 0x00F0)>>4;

	const float x1min = pha_range[0][0];
	const float x2min = pha_range[1][0];

	const float rdx1 = pha_nx[0] / ( pha_range[0][1] - pha_range[0][0] );
	const float rdx2 = pha_nx[1] / ( pha_range[1][1] - pha_range[1][0] );

	for ( int i = 0; i<spec->np; i+=BUF_SIZE ) {
		int np = ( i + BUF_SIZE > spec->np )? spec->np - i : BUF_SIZE;

		spec_pha_axis( spec, i, np, quant1, pha_x1 );
	    spec_pha_axis( spec, i, np, quant2, pha_x2 );

		for ( int k = 0; k < np; k++ ) {

			float nx1 = ( pha_x1[k] - x1min ) * rdx1;
			float nx2 = ( pha_x2[k] - x2min ) * rdx2;

			int i1 = (int)(nx1 + 0.5f);
			int i2 = (int)(nx2 + 0.5f);

			float w1 = nx1 - i1 + 0.5f;
			float w2 = nx2 - i2 + 0.5f;

			int idx = i1 + nrow*i2;

			if ( i2 >= 0 && i2 < pha_nx[1] ) {

				if (i1 >= 0 && i1 < pha_nx[0]) {
					buf[ idx ] += (1.0f-w1)*(1.0f-w2)*spec->q;
				}

				if (i1+1 >= 0 && i1+1 < pha_nx[0] ) {
					buf[ idx + 1 ] += w1*(1.0f-w2)*spec->q;
				}
			}

			idx += nrow;
			if ( i2+1 >= 0 && i2+1 < pha_nx[1] ) {

				if (i1 >= 0 && i1 < pha_nx[0]) {
					buf[ idx ] += (1.0f-w1)*w2*spec->q;
				}

				if (i1+1 >= 0 && i1+1 < pha_nx[0] ) {
					buf[ idx + 1 ] += w1*w2*spec->q;
				}
			}

		}

	}
}

void spec_rep_pha( const t_species *spec, const int rep_type,
			  const int pha_nx[], const float pha_range[][2] )
{

	char const * const pha_ax_name[] = {"x1","x2","x3","u1","u2","u3"};
	char pha_name[64];

	// Allocate phasespace buffer
	float* restrict buf = malloc( pha_nx[0] * pha_nx[1] * sizeof( float ));
	memset( buf, 0, pha_nx[0] * pha_nx[1] * sizeof( float ));

	// Deposit the phasespace
	spec_deposit_pha( spec, rep_type, pha_nx, pha_range, buf );

	// save the data in hdf5 format
	int quant1 = rep_type & 0x000F;
	int quant2 = (rep_type & 0x00F0)>>4;

    const char * pha_ax1_units = spec_pha_axis_units(quant1);
    const char * pha_ax2_units = spec_pha_axis_units(quant2);

	sprintf( pha_name, "%s%s", pha_ax_name[quant1-1], pha_ax_name[quant2-1] );

    t_zdf_grid_axis axis[2];
    axis[0] = (t_zdf_grid_axis) {
    	.min = pha_range[0][0],
    	.max = pha_range[0][1],
    	.label = (char *) pha_ax_name[ quant1 - 1 ],
    	.units = (char *) pha_ax1_units
    };

    axis[1] = (t_zdf_grid_axis) {
    	.min = pha_range[1][0],
    	.max = pha_range[1][1],
    	.label = (char *) pha_ax_name[ quant2 - 1 ],
    	.units = (char *) pha_ax2_units
    };

    t_zdf_grid_info info = {
    	.ndims = 2,
    	.label = pha_name,
    	.units = "a.u.",
    	.axis  = axis
    };

    info.nx[0] = pha_nx[0];
    info.nx[1] = pha_nx[1];

    t_zdf_iteration iter = {
    	.n = spec->iter,
    	.t = spec -> iter * spec -> dt,
    	.time_units = "1/\\omega_p"
    };

	zdf_save_grid( buf, &info, &iter, spec->name );

	// Free temp. buffer
	free( buf );

}


void spec_report( const t_species *spec, const int rep_type, 
				  const int pha_nx[], const float pha_range[][2] )
{
	
	switch (rep_type & 0xF000) {
		case CHARGE:
			spec_rep_charge( spec );
			break;

		case PHA:
			spec_rep_pha( spec, rep_type, pha_nx, pha_range );
			break;

		case PARTICLES:
			spec_rep_particles( spec );
			break;
	}
	
	
}
