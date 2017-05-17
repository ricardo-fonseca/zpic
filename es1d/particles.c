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

#include "grid.h"

#include "particles.h"

#include "random.h"
#include "field.h"
#include "charge.h"

#include "zdf.h"
#include "timer.h"

static double _spec_time = 0.0;

void spec_sort( t_species *spec );

/*********************************************************************************************

 Initialization

 *********************************************************************************************/


double spec_time()
{
	return _spec_time;
}


void spec_set_u( t_species* spec, const int start, const int end )
{
	int i;

	for (i = start; i <= end; i++) {
		spec->part[i].vx = spec -> vfl + spec -> vth * rand_norm();
	}

}

void spec_set_x( t_species* spec, const int range[] )
{

	int i, k, ip;

	float edge, start, end;

	// Calculate particle positions inside the cell
	const int npc = spec->ppc;

	float poscell[npc];

	for (i=0; i<spec->ppc; i++) {
		poscell[i]   = (1 + 2*i - npc) / (2.0*npc);
	}

	ip = spec -> np;

	// Set position of particles in the specified grid range according to the density profile
	switch ( spec -> density.type ) {
	case STEP: // Step like density profile

		// Get edge position normalized to cell size;
		edge = spec -> density.edge / spec -> dx - 0.5;

		for (i = range[0]; i <= range[1]; i++) {

			for (k=0; k<npc; k++) {
				if ( i + poscell[k] > edge ) {
					spec->part[ip].ix = i;
					spec->part[ip].x = poscell[k];
					ip++;
				}
			}
		}
		break;

	case SLAB: // Step like density profile

		// Get edge position normalized to cell size;
		start = spec -> density.start / spec -> dx - 0.5;
		end   = spec -> density.end / spec -> dx - 0.5;

		for (i = range[0]; i <= range[1]; i++) {

			for (k=0; k<npc; k++) {
				if ( i + poscell[k] > start &&  i + poscell[k] < end ) {
					spec->part[ip].ix = i;
					spec->part[ip].x = poscell[k];
					ip++;
				}
			}
		}
		break;

	default: // Uniform density
		for (i = range[0]; i <= range[1]; i++) {

			for (k=0; k<npc; k++) {
				spec->part[ip].ix = i;
				spec->part[ip].x = poscell[k];
				ip++;
			}
		}
	}

	spec -> np = ip;

}

void spec_inject_particles( t_species* spec, const int range[] )
{
	int start = spec -> np;

	// Get maximum number of particles to inject
	int np_inj = ( range[1] - range[0] + 1 ) * spec -> ppc;

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

void spec_new( t_species* spec, char name[], const float m_q, const int ppc,
			  const float *vfl, const float * vth,
			  const int nx, float box, const float dt, t_density* density )
{

	int npc;

	// Species name
	strncpy( spec -> name, name, MAX_SPNAME_LEN );


	spec->nx = nx;
	spec->ppc = ppc;
	npc = ppc;

	spec->box = box;
	spec->dx = box / nx;

	spec -> m_q = m_q;
	spec -> q = copysign( 1.0f, m_q ) / npc;

	spec -> dt = dt;

	// Initialize particle buffer
	spec->np_max = 0;
	spec->part = NULL;


	// Initialize density profile
	if ( density ) {
		spec -> density = *density;
	} else {
		spec -> density.type = UNIFORM;
	}

	// Initialize temperature profile
	if ( vfl ) {
		spec -> vfl = *vfl;
	} else {
		spec -> vfl = 0;
	}

	if ( vth ) {
		spec -> vth = *vth;
	} else {
		spec -> vth = 0;
	}

	// Reset iteration number
	spec -> iter = 0;

    // Inject initial particle distribution
    spec -> np = 0;

    const int range[2] = {0, nx-1};

    spec_inject_particles( spec, range );

}

void spec_delete( t_species* spec )
{
	free(spec->part);
	spec->np = -1;
}


/*********************************************************************************************

 Charge deposition

 *********************************************************************************************/

void deposit( t_scalar_grid * rho, const t_part* restrict const part, float q )
{
	int i;
	float s0, s1;

	i = part->ix;

	s0 = 0.5f - part->x;
	s1 = 0.5f + part->x;

	rho->s[i] += s0 * q;
	rho->s[i+1] += s1 * q;

}

/*********************************************************************************************

 Sorting

 *********************************************************************************************/

void spec_sort( t_species* spec )
{
	int *idx, *npic;

	int ncell = spec->nx;

	// Allocate index memory
	idx  = malloc(spec->np*sizeof(int));

	// Allocate temp. array with number of particles in cell
	npic = malloc( ncell * sizeof(int));
	memset( npic, 0, ncell * sizeof(int));

	// Generate sorted index
    int i;
	for (i=0; i<spec->np; i++) {
		idx[i] = spec->part[i].ix;
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


void interpolate_fld( t_scalar_grid* E,
	          const t_part* restrict const part, float* restrict const Ex )
{
	register int i;
	register float s0, s1;

	i = part->ix;

	s0 = 0.5f - part->x;
	s1 = 0.5f + part->x;

	*Ex = E->s[i] * s0 + E->s[i+1] * s1;
}

int ltrim( float x )
{
	// If the particles do not move more than one cell we can use this
	// return (( x >= 0.5f )?1:0) - (( x < -0.5f )?1:0);

	// To allow motion over more than 1 cell we must change the function
	int di = x + ((x>0)?0.5f:-0.5f);

	// Correct for values lying on the left cell boundary
	di += (x-di == 0.5f);

	return di;
}

void spec_advance( t_species* spec, t_field* field, t_charge* charge )
{
	int i;

	uint64_t t0;
	t0 = timer_ticks();

	const float tem   = spec->dt / spec -> m_q;
	const float dt_dx = spec->dt / spec->dx;

	const int nx0 = spec -> nx;

	double energy = 0;

	// Advance particles
	for (i=0; i<spec->np; i++) {

		float Ex, vx;


		float x1;

		int di;
		float dx;

		// Load particle velocity
		vx = spec -> part[i].vx;

		// interpolate field
		interpolate_fld( &field -> E, &spec -> part[i], &Ex );

		// advance u / x using a leap frog scheme
		vx += tem * Ex;

		// Accumulate particle energy K = 1/2 (vold * vnew)
		// (the 1/2 was moved out of the loop)
		energy += spec -> part[i].vx * vx;

		// Store new momenta
		spec -> part[i].vx = vx;

		// push particle
		dx = dt_dx * vx;

		x1 = spec -> part[i].x + dx;

		di = ltrim(x1);

		x1 -= di;

		// Store results
		spec -> part[i].x = x1;
		spec -> part[i].ix += di;

		// Use periodic boundaries in x
		spec -> part[i].ix += (( spec -> part[i].ix < 0    ) ? nx0 : 0 ) -
		                      (( spec -> part[i].ix >= nx0 ) ? nx0 : 0 );

		// Deposit charge for the particle (t+1)
		deposit( &charge -> rho, &spec -> part[i], spec->q );

	}

	// Store energy
	spec -> energy = 0.5 * spec-> q * spec -> m_q * energy * spec -> dx;

	// Advance internal iteration number
  spec -> iter += 1;
	
	// Sort species at every 16 time steps
	if ( ! (spec -> iter % 16) ) spec_sort( spec );

	_spec_time += timer_interval_seconds( t0, timer_ticks() );
}


/*********************************************************************************************

 Charge Deposition

 *********************************************************************************************/


void spec_deposit_charge( const t_species* spec, float* charge )
{
	int i;

	// Charge array is expected to have 1 guard cell at the upper boundary

	float q = spec -> q;

	for (i=0; i<spec->np; i++) {
		int idx = spec->part[i].ix;
		float w1 = spec->part[i].x;

		charge[ idx            ] += ( 0.5f - w1 ) * q;
		charge[ idx + 1        ] += ( 0.5f + w1 ) * q;
	}

}

/*********************************************************************************************

 Diagnostics

 *********************************************************************************************/

void spec_rep_particles( const t_species *spec )
{

	t_zdf_file part_file;

	unsigned i;

	const char * quants[] = {
	    "x1",
	    "v1"
	};

	const char * units[] = {
	    "c/\\omega_p",
	    "c"
	};

    t_zdf_iteration iter = {
    	.n = spec->iter,
    	.t = spec -> iter * spec -> dt,
    	.time_units = "1/\\omega_p"
    };

	// Allocate buffer for positions

	t_zdf_part_info info = {
		.name = (char *) spec -> name,
		.nquants = 2,
		.quants = (char **) quants,
		.units = (char **) units,
		.np = spec ->np
	};

	// Create file and add description
	zdf_part_file_open( &part_file, &info, &iter, "PARTICLES" );

	// Add positions and generalized velocities
	size_t size = ( spec -> np ) * sizeof( float );
	float* data = malloc( size );

	// x
	for( i = 0; i < spec ->np; i++ )
		data[i] = (spec -> part[i].ix + (spec -> part[i].x + 0.5f)) * spec -> dx;
	zdf_part_file_add_quant( &part_file, quants[0], data, spec ->np );

	// vx
	for( i = 0; i < spec ->np; i++ ) data[i] = spec -> part[i].vx;
	zdf_part_file_add_quant( &part_file, quants[1], data, spec ->np );

	free( data );

	zdf_close_file( &part_file );
}


void spec_rep_charge( const t_species *spec )
{
	float *charge;
	size_t size;
	int i;

	// Add 1 guard cell to the upper boundary
	size = ( spec -> nx + 1 ) * sizeof( float );
	charge = malloc( size );
	memset( charge, 0, size );

	// Deposit the charge
	spec_deposit_charge( spec, charge );

	// Correct boundary values - x
	charge[ 0 ] += charge[ spec -> nx ];

	// Compact the data to save the file (throw away guard cells)
	float buffer[ spec -> nx ];

	for ( i = 0; i < spec->nx; i++ ) {
		buffer[i] = charge[i];
	}

	free( charge );

    t_zdf_grid_axis axis = {
    	.min = 0.0,
    	.max = spec->box,
    	.label = "x_1",
    	.units = "c/\\omega_p"
    };

    t_zdf_grid_info info = {
    	.ndims = 1,
    	.label = "charge",
    	.units = "n_e",
    	.axis  = &axis
    };

    info.nx[0] = spec->nx;

    t_zdf_iteration iter = {
    	.n = spec->iter,
    	.t = spec -> iter * spec -> dt,
    	.time_units = "1/\\omega_p"
    };

	zdf_save_grid( buffer, &info, &iter, spec->name );

}


void spec_pha_axis( const t_species *spec, int i0, int np, int quant, float *axis )
{
	int i;

	switch (quant) {
		case X1:
			for (i = 0; i < np; i++)
				axis[i] = ( (spec -> part[i0+i].x + 0.5f) + spec -> part[i0+i].ix ) * spec -> dx;
			break;
		case V1:
			for (i = 0; i < np; i++)
				axis[i] = spec -> part[i0+i].vx;
			break;
	}
}

const char * spec_pha_axis_units( int quant ) {
	switch (quant) {
		case X1:
			return("c/\\omega_p");
			break;
		case V1:
			return("c");
	}
	return("");
}


void spec_rep_pha( const t_species *spec, const int rep_type,
			  const int pha_nx[], const float pha_range[][2] )
{
	const int BUF_SIZE = 1024;
	float pha_x1[BUF_SIZE], pha_x2[BUF_SIZE];

	int i, nrow;

	int quant1, quant2;
	float rdx1, rdx2, x1min, x2min;

	char const * const pha_ax_name[] = {"x1","x2","x3","v1","v2","v3"};
	char pha_name[64];

	// Allocate phasespace buffer
	float* restrict buf = malloc( pha_nx[0] * pha_nx[1] * sizeof( float ));
	memset( buf, 0, pha_nx[0] * pha_nx[1] * sizeof( float ));

	nrow = pha_nx[0];

	quant1 = rep_type & 0x000F;
	quant2 = (rep_type & 0x00F0)>>4;

    const char * pha_ax1_units = spec_pha_axis_units(quant1);
    const char * pha_ax2_units = spec_pha_axis_units(quant2);

	x1min = pha_range[0][0];
	x2min = pha_range[1][0];

	rdx1 = pha_nx[0] / ( pha_range[0][1] - pha_range[0][0] );
	rdx2 = pha_nx[1] / ( pha_range[1][1] - pha_range[1][0] );

	for (i = 0; i<spec->np; i+=BUF_SIZE) {
		int k;
		int np = ( i + BUF_SIZE > spec->np )? spec->np - i : BUF_SIZE;

		spec_pha_axis( spec, i, np, quant1, pha_x1 );
	    spec_pha_axis( spec, i, np, quant2, pha_x2 );

		for (k = 0; k < np; k++) {
			float nx1, nx2, w1, w2;
			int i1, i2, idx;

			nx1 = ( pha_x1[k] - x1min ) * rdx1;
			nx2 = ( pha_x2[k] - x2min ) * rdx2;

			i1 = (int)(nx1 + 0.5f);
			i2 = (int)(nx2 + 0.5f);

			w1 = nx1 - i1 + 0.5f;
			w2 = nx2 - i2 + 0.5f;

			idx = i1 + nrow*i2;

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

	// save the data in hdf5 format
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
