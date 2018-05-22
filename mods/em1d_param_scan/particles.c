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

#include "particles.h"

#include "random.h"
#include "emf.h"
#include "current.h"

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
		spec->part[i].ux = spec -> ufl[0] + spec -> uth[0] * rand_norm();
		spec->part[i].uy = spec -> ufl[1] + spec -> uth[1] * rand_norm();
		spec->part[i].uz = spec -> ufl[2] + spec -> uth[2] * rand_norm();
	}

}

/**
 * Number of particles to be injected.
 *
 * Calculates the number of particles to be injected in the specified range according
 * to the specified density profile. The returned value is not exact but it is
 * guaranteed to be larger than the actual number of particles to be injected
 *
 * @param spec Particle species
 * @param range[] Range of cells in which to inject
 * @return Number of particles to be injected
 */
int spec_np_inj( t_species* spec, const int range[] )
{
	int np_inj;

	switch ( spec -> density.type ) {
	case STEP: // Step like density profile
		{
			int i0 = spec -> density.start / spec -> dx - spec -> n_move;

			if ( i0 > range[1] ) {
				np_inj = 0;
			} else {
				if ( i0 < range[0] ) i0 = range[0];
				np_inj = ( range[1] - i0 + 1 ) * spec -> ppc;
			}
		}
		break;

	case SLAB: // Slab like density profile
		{
			int i0 = spec -> density.start / spec -> dx - spec -> n_move;
			int i1 = spec -> density.end / spec -> dx - spec -> n_move;

			if ( (i0 > range[1]) || (i1 < range[0]) ) {
				np_inj = 0;
			} else {
				if ( i0 < range[0] ) i0 = range[0];
				if ( i1 > range[1] ) i1 = range[1];
				np_inj = ( i1 - i0 + 1 ) * spec -> ppc;
			}
		}
		break;

	case RAMP: // ramp density profile
		{
			// Ramp start / finish
			float x0 = spec -> density.start;
			float x1 = spec -> density.end;

			// Injection start / finish positions (in simulation units)
			float a = (range[0] + spec -> n_move) * spec->dx;
			float b = (range[1] + 1 + spec -> n_move) * spec->dx;

			// If outside of ramp (or invalid ramp) return 0
			if ( (x1 <= x0) || (a > x1) || (b < x0) ) {
				np_inj = 0;
			} else {
				// limit integration boundaries to ramp start/end
				if ( a < x0 ) a = x0;
				if ( b > x1 ) b = x1;

				// Get total injected charge
				float n0 = spec -> density.ramp[0];
				float n1 = spec -> density.ramp[1];
				float q = (b-a)*( n0 + 0.5 * (a+b-2*x0)*(n1-n0)/(x1-x0));

				// Get corresponding number of simulation particles
				np_inj = q * spec -> ppc / spec -> dx;

			}
		}
		break;

	case CUSTOM: // custom density profile
		{

			// Integrate total charge
			double q = 0.5 * ( (*spec -> density.custom)((range[0] + spec -> n_move) * spec->dx) +
				               (*spec -> density.custom)((range[1] + spec -> n_move) * spec->dx) );

			for( int i = range[0]+1; i < range[1]; i++) {
				q += (*spec -> density.custom)((i + spec -> n_move) * spec->dx);
			}

			// Get corresponding number of simulation particles, rounding up
			np_inj = ceil(q * spec -> ppc);

		}
		break;

	default: // Uniform density
		np_inj = ( range[1] - range[0] + 1 ) * spec -> ppc;
	}

	// printf("Predicts injecting %d particles\n", np_inj);
	return np_inj;

}

void spec_set_x( t_species* spec, const int range[] )
{

	int i, k, ip;

	float start, end;

	// Calculate particle positions inside the cell
	const int npc = spec->ppc;

	float poscell[npc];

	for (i=0; i<spec->ppc; i++) {
		poscell[i]   = ( i + 0.5 ) / npc;
	}

	ip = spec -> np;

	// Set position of particles in the specified grid range according to the density profile
	switch ( spec -> density.type ) {
	case STEP: // Step like density profile

		// Get edge position normalized to cell size;
		start = spec -> density.start / spec -> dx - spec -> n_move;

		for (i = range[0]; i <= range[1]; i++) {

			for (k=0; k<npc; k++) {
				if ( i + poscell[k] > start ) {
					spec->part[ip].ix = i;
					spec->part[ip].x = poscell[k];
					ip++;
				}
			}
		}

		// printf("Injected %d particles with step injection \n", ip - spec -> np );
		break;

	case SLAB: // Slab like density profile

		// Get edge position normalized to cell size;
		start = spec -> density.start / spec -> dx - spec -> n_move;
		end   = spec -> density.end / spec -> dx - spec -> n_move;

		for (i = range[0]; i <= range[1]; i++) {

			for (k=0; k<npc; k++) {
				if ( i + poscell[k] > start &&  i + poscell[k] < end ) {
					spec->part[ip].ix = i;
					spec->part[ip].x = poscell[k];
					ip++;
				}
			}
		}

		// printf("Injected %d particles with slab injection \n", ip - spec -> np );
		break;

	case RAMP: // ramp like density profile

		{
			// Ramp start/finish in cell units
			double r0 = spec -> density.start / spec -> dx;
			double r1 = spec -> density.end / spec -> dx;

			// If outside ramp return
			if (((range[0] + spec -> n_move) > r1 ) ||
				((range[1] + spec -> n_move) < r0 )) break;


			double n0 = spec -> density.ramp[0];
			double n1 = spec -> density.ramp[1];

			// Only consider the ramp for x > 0
			if ( r0 < 0 ) {
				n0 += - r0 * (n1-n0) / (r1-r0);
				r0 = 0;
			}

            // Charge per simulation particle
			double cpp = 1.0 / spec->ppc;

			for( k = spec -> density.total_np_inj; ; k++ ) {
				// Desired cumulative density, normalized to the [0,1] interval
				double Rs = (k+0.5) * cpp / (r1 - r0);

				// Position normalized to the [0,1] interval
				// double pos = (-a + sqrt( a*a + 2 * b * Rs ))/b;
				double pos = 2 * Rs / (sqrt( n0*n0 + 2 * (n1-n0) * Rs ) + n0);

				// If outside of ramp interval we are done
				if ( pos > 1 ) break;

				// Position in simulation cell units
				pos = r0 + (r1-r0) * pos;

				// Injection cell
				int ix = pos;

				// (*debug*) This must never happen
				if ( ix - spec -> n_move < range[0] ) {
					fprintf(stderr, "(*error*) attempting to inject outside of valid range.\n");
					break;
				}

				// If outside injection range we are done
				if ( ix - spec -> n_move > range[1] ) break;

				// Inject particle
				spec->part[ip].ix = ix - spec -> n_move;
				spec->part[ip].x = pos - ix;
				ip++;

			}
        }

		// printf("Injected %d particles with ramp injection \n", ip - spec -> np );
		break;

	case CUSTOM: // ramp like density profile

		{

			const double dx = spec -> dx;

			// Charge per simulation particle
			const double cpp = 1.0 / spec->ppc;

			// Injected particles
			k = spec -> density.total_np_inj;

			int ix = range[0];

			// Density on cell edges
			double n0;
			double n1 = (*spec -> density.custom)((ix + spec -> n_move) * dx);

			// Accumulated density on cell edges
			double d0;
			double d1 = spec -> density.custom_q_inj;

			double Rs;

			while( ix <= range[1] ){

				// Get density on the edges of current cell
				n0 = n1;
				n1 = (*spec -> density.custom)((ix + 1 + spec -> n_move)*dx);

				// Get cumulative density on the edges of current cell
				d0 = d1;
				d1 += 0.5 * (n0+n1);

				while( ( Rs =  (k+0.5) * cpp ) < d1 ) {

					// Quadratic formula
					// double pos = (-n0 + sqrt( n0*n0 + 2 * (n1-n0) * (Rs-d0) ))/(n1-n0);

					// This version avoids a division by 0 if n1 = n0
					double pos = 2 * (Rs-d0) /( sqrt( n0*n0 + 2 * (n1-n0) * (Rs-d0) ) + n0 );

					spec->part[ip].ix = ix;
					spec->part[ip].x = pos;
					ip++;

					k++;
				}

				// Move to next cell
				ix++;
			}

			spec -> density.custom_q_inj = d1;
        }

		// printf("Injected %d particles with custom injection \n", ip - spec -> np );
		break;

	default: // Uniform density
		for (i = range[0]; i <= range[1]; i++) {

			for (k=0; k<npc; k++) {
				spec->part[ip].ix = i;
				spec->part[ip].x = poscell[k];
				ip++;
			}
		}
		// printf("Injected %d particles with uniform injection \n", ip - spec -> np );

	}

	// Update total number of injected particles
	spec -> density.total_np_inj += ip - spec -> np;

	// Update number of particles in buffer
	spec -> np = ip;

}

void spec_inject_particles( t_species* spec, const int range[] )
{
	int start = spec -> np;

	// Get maximum number of particles to inject
	int np_inj = spec_np_inj( spec, range );

	// Check if buffer is large enough and if not reallocate
	if ( spec -> np + np_inj > spec -> np_max ) {
        // Increase by chunks of 1024 particles
        spec -> np_max = (( spec -> np_max + np_inj )/1024 + 1) * 1024;
		spec -> part = realloc( (void*) spec -> part, spec -> np_max * sizeof(t_part) );
	}

	// Set particle positions
	spec_set_x( spec, range );

	// Set momentum of injected particles
	spec_set_u( spec, start, spec -> np - 1 );

}

void spec_new( t_species* spec, char name[], const float m_q, const int ppc,
			  const float *ufl, const float * uth,
			  const int nx, float box, const float dt, t_density* density )
{

	int i, npc;

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
		if ( spec -> density.n == 0. ) spec -> density.n = 1.0;
	} else {
		// Default values
		spec -> density = (t_density) { .type = UNIFORM, .n = 1.0 };
	}
	spec -> density.total_np_inj = 0;
	spec -> density.custom_q_inj = 0.;

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

	// Reset moving window information
	spec -> moving_window = 0;
	spec -> n_move = 0;

    // Inject initial particle distribution
    spec -> np = 0;

    const int range[2] = {0, nx-1};

    spec_inject_particles( spec, range );

}

void spec_move_window( t_species *spec ){

	if ((spec->iter * spec->dt ) > (spec->dx * (spec->n_move + 1)))  {

        // shift all particles left
        // particles leaving the box will be removed later
        int i;
        for( i = 0; i < spec->np; i++ ) {
        	spec->part[i].ix--;
        }

		// Increase moving window counter
		spec -> n_move++;

        // Inject particles in the right edge of the simulation box
	    const int range[2] = {spec->nx-1,spec->nx-1};
	    spec_inject_particles( spec, range );

	}

}

void spec_delete( t_species* spec )
{
	free(spec->part);
	spec->np = -1;
}


/*********************************************************************************************

 Cuurent deposition

 *********************************************************************************************/

void dep_current_esk( int ix0, int di,
						float x0, float x1,
						float qnx, float qvy, float qvz,
						t_current *current )
{

	int i;
	t_fld S0x[4], S1x[4], DSx[4];
	t_fld Wx[4], Wy[4], Wz[4];

	S0x[0] = 0.0f;
	S0x[1] = 1.0f - x0;
	S0x[2] = x0;
	S0x[3] = 0.0f;

	for (i=0; i<4; i++) {
		S1x[i] = 0.0f;
	}

	S1x[ 1 + di ] = 1.0f - x1;
	S1x[ 2 + di ] = x1;

	for (i=0; i<4; i++) {
		DSx[i] = S1x[i] - S0x[i];
	}

	for (i=0; i<4; i++) {
		Wx[i] = qnx * DSx[i];
		Wy[i] = qvy * (S0x[i] + DSx[i]/2.0f);
		Wz[i] = qvz * (S0x[i] + DSx[i]/2.0f);
	}

	t_vfld* restrict const J = current -> J;
	// jx
	t_fld c;

	c = - Wx[0];
	J[ ix0 - 1 ].x += c;
	for (i=1; i<4; i++) {
		c -=  Wx[i];
		J[ ix0 + i ].x += c;
	}

	// jy, jz
	for (i=0; i<4; i++) {
		J[ ix0 + i - 1 ].y += Wy[ i ];
		J[ ix0 + i - 1 ].z += Wz[ i ];
	}


}

void dep_current_zamb( int ix0, int di,
						float x0, float dx,
						float qnx, float qvy, float qvz,
						t_current *current )
{
	// Split the particle trajectory

	typedef struct {
		float x0, x1, dx, qvy, qvz;
		int ix;
	} t_vp;

	t_vp vp[3];
	int vnp = 1;

	// split
	vp[0].x0 = x0;
	vp[0].dx = dx;

	vp[0].x1 = x0+dx;

	vp[0].qvy = qvy/2.0;
	vp[0].qvz = qvz/2.0;

	vp[0].ix = ix0;

	// x split
	if ( di != 0 ) {

		//int ib = ( di+1 )>>1;
		int ib = ( di == 1 );

		float delta = (x0+dx-ib)/dx;

		// Add new particle
		vp[1].x0 = 1-ib;
		vp[1].x1 = (x0 + dx) - di;
		vp[1].dx = dx*delta;
		vp[1].ix = ix0 + di;


		vp[1].qvy = vp[0].qvy*delta;
		vp[1].qvz = vp[0].qvz*delta;

		// Correct previous particle
		vp[0].x1 = ib;
		vp[0].dx *= (1.0f-delta);

		vp[0].qvy *= (1.0f-delta);
		vp[0].qvz *= (1.0f-delta);

        vnp++;
	}

	// Deposit virtual particle currents
	int k;
	t_vfld* restrict const J = current -> J;

	for (k = 0; k < vnp; k++) {
		float S0x[2], S1x[2];

		S0x[0] = 1.0f - vp[k].x0;
		S0x[1] = vp[k].x0;

		S1x[0] = 1.0f - vp[k].x1;
		S1x[1] = vp[k].x1;

		J[ vp[k].ix ].x += qnx * vp[k].dx;

		J[ vp[k].ix     ].y += vp[k].qvy * (S0x[0]+S1x[0]+(S0x[0]-S1x[0])/2.0f);
		J[ vp[k].ix + 1 ].y += vp[k].qvy * (S0x[1]+S1x[1]+(S0x[1]-S1x[1])/2.0f);

		J[ vp[k].ix     ].z += vp[k].qvz * (S0x[0]+S1x[0]+(S0x[0]-S1x[0])/2.0f);
		J[ vp[k].ix  +1 ].z += vp[k].qvz * (S0x[1]+S1x[1]+(S0x[1]-S1x[1])/2.0f);
	}

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
		t_part tmp;
		int k;

		k = idx[i];
		while ( k > i ) {
			int t;

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


void interpolate_fld( const t_vfld* restrict const E, const t_vfld* restrict const B,
	          const t_part* restrict const part, t_vfld* restrict const Ep, t_vfld* restrict const Bp )
{
	int i, ih;
	t_fld w1, w1h;

	i = part->ix;

	w1 = part->x;
	ih = (w1 <0.5f)? -1 : 0;
	w1h = w1 + ((w1 <0.5f)?0.5f:-0.5f);


	ih += i;

	Ep->x = E[ih].x * (1.0f - w1h) + E[ih+1].x * w1h;
	Ep->y = E[i ].y * (1.0f -  w1) + E[i+1 ].y * w1;
	Ep->z = E[i ].z * (1.0f -  w1) + E[i+1 ].z * w1;

	Bp->x = B[i ].x * (1.0f  - w1) + B[i+1 ].x * w1;
	Bp->y = B[ih].y * (1.0f - w1h) + B[ih+1].y * w1h;
	Bp->z = B[ih].z * (1.0f - w1h) + B[ih+1].z * w1h;

}

int ltrim( float x )
{
	return ( x >= 1.0f ) - ( x < 0.0f );
}

void spec_advance( t_species* spec, t_emf* emf, t_current* current )
{
	int i;
	float qnx, qvy, qvz;

	uint64_t t0;
	t0 = timer_ticks();

	const float tem   = 0.5 * spec->dt/spec -> m_q;
	const float dt_dx = spec->dt / spec->dx;

	// Auxiliary values for current deposition
	qnx = spec -> q *  spec->dx / spec->dt;

	const int nx0 = spec -> nx;

	double energy = 0;

	// Advance particles
	for (i=0; i<spec->np; i++) {

		t_vfld Ep, Bp;
		float utx, uty, utz;
		float ux, uy, uz, u2;
		float gamma, rg, gtem, otsq;

		float x1;

		int di;
		float dx;

		// Load particle momenta
		ux = spec -> part[i].ux;
		uy = spec -> part[i].uy;
		uz = spec -> part[i].uz;

		// interpolate fields
		interpolate_fld( emf -> E, emf -> B, &spec -> part[i], &Ep, &Bp );
		// Ep.x = Ep.y = Ep.z = Bp.x = Bp.y = Bp.z = 0;

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

		x1 = spec -> part[i].x + dx;

		di = ltrim(x1);

		x1 -= di;

		qvy = spec->q * uy * rg;
		qvz = spec->q * uz * rg;

		// deposit current using Eskirepov method
		// dep_current_esk( spec -> part[i].ix, di,
		// 				 spec -> part[i].x, x1,
		// 				 qnx, qvy, qvz,
		// 				 current );

		dep_current_zamb( spec -> part[i].ix, di,
						 spec -> part[i].x, dx,
						 qnx, qvy, qvz,
						 current );

		// Store results
		spec -> part[i].x = x1;
		spec -> part[i].ix += di;

	}

	// Advance internal iteration number
    spec -> iter += 1;

    // Check for particles leaving the box
	if ( spec -> moving_window ){

        // Move simulation window if needed
		spec_move_window( spec );

		// Use absorbing boundaries along x
		i = 0;
		while ( i < spec -> np ) {
			if (( spec -> part[i].ix < 0 ) || ( spec -> part[i].ix >= nx0 )) {
				spec -> part[i] = spec -> part[ -- spec -> np ];
				continue;
			}
			i++;
		}

	} else {
		// Use periodic boundaries in x
		for (i=0; i<spec->np; i++) {
			spec -> part[i].ix += (( spec -> part[i].ix < 0 ) ? nx0 : 0 ) - (( spec -> part[i].ix >= nx0 ) ? nx0 : 0);
		}
	}

	// Sort species at every 16 time steps
	if ( ! (spec -> iter % 16) ) spec_sort( spec );


    // Move simulation window if needed
    if ( spec -> moving_window )


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
		t_fld w1 = spec->part[i].x;

		charge[ idx            ] += ( 1.0f - w1 ) * q;
		charge[ idx + 1        ] += (        w1 ) * q;
	}

	// Correct boundary values

	// x
	if ( ! spec -> moving_window ){
		charge[ 0 ] += charge[ spec -> nx ];
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
	    "u1","u2","u3"
	};

	const char * units[] = {
	    "c/\\omega_p",
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
		.nquants = 4,
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
		data[i] = (spec -> n_move + spec -> part[i].ix + spec -> part[i].x ) * spec -> dx;
	zdf_part_file_add_quant( &part_file, quants[0], data, spec ->np );

	// ux
	for( i = 0; i < spec ->np; i++ ) data[i] = spec -> part[i].ux;
	zdf_part_file_add_quant( &part_file, quants[1], data, spec ->np );

	// uy
	for( i = 0; i < spec ->np; i++ ) data[i] = spec -> part[i].uy;
	zdf_part_file_add_quant( &part_file, quants[2], data, spec ->np );

	// uz
	for( i = 0; i < spec ->np; i++ ) data[i] = spec -> part[i].uz;
	zdf_part_file_add_quant( &part_file, quants[3], data, spec ->np );

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

	// Compact the data to save the file (throw away guard cells)
	float buffer[ spec -> nx ];

	for ( i = 0; i < spec->nx; i++ ) {
		buffer[i] = charge[i];
	}

	free( charge );

	// Set grid boundaries accounting for moving window
    t_zdf_grid_axis axis = {
    	.min = spec -> n_move * spec -> dx,
    	.max = spec->box + spec -> n_move * spec -> dx,
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
				axis[i] = ( spec -> part[i0+i].x + spec -> part[i0+i].ix ) * spec -> dx;
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


void spec_rep_pha( const t_species *spec, const int rep_type,
			  const int pha_nx[], const float pha_range[][2] )
{
	const int BUF_SIZE = 1024;
	float pha_x1[BUF_SIZE], pha_x2[BUF_SIZE];

	int i, nrow;

	int quant1, quant2;
	float rdx1, rdx2, x1min, x2min;

	char const * const pha_ax_name[] = {"x1","x2","x3","u1","u2","u3"};
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

