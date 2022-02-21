/**
 * ZPIC - em2d
 *
 * External field initialization options
 */

#include <stdlib.h>
#include <math.h>
#include "../simulation.h"

/**
 * Field from an infinite wire with unitary current
 * perpendicular to the simulation plane
 */


float3 ext_B( int ix, float dx, int iy, float dy, void *data ) {
    // Wire position
    const float x0 = 6.4;
    const float y0 = 6.4;

	float3 b;
    float x, y, r2;

    // Field components must be calculated in the correct positions
    // inside the cells

    // Bx
    x = ix*dx       - x0;
    y = (iy+0.5)*dy - y0;

    r2 = x*x+y*y;
    b.x = -y/r2;

    // By
    x = (ix+0.5)*dx - x0;
    y = iy*dy       - y0;

    r2 = x*x+y*y;
    b.y = x/r2;

    // Bz
    // x = (ix+0.5)*dx - x0;
    // y = (ix+0.5)*dy - y0;
    b.z = 0;

    return b;
}

void sim_init( t_simulation* sim ){

	// Time step
	float dt = 0.07;
	float tmax = 0.14;

	// Simulation box
	int   nx[2]  = { 128, 128 };
	float box[2] = { 12.8, 12.8 };

	// Diagnostic frequency
	int ndump = 1;

    // Initialize particles
	const int n_species = 1;

	// Use 2x2 particles per cell
	int ppc[] = {8,8};

	t_species* species = (t_species *) malloc( n_species * sizeof( t_species ));

	spec_new( &species[0], "electrons", -1.0, ppc, NULL, NULL, nx, box, dt, NULL );

	// Initialize Simulation data
	sim_new( sim, nx, box, dt, tmax, ndump, species, n_species );

    t_emf_ext_fld ext_fld = {
		.B_type = EMF_FLD_TYPE_CUSTOM,
		.B_custom = &ext_B,
	};
	sim_set_ext_fld( sim, &ext_fld );

}

void sim_report( t_simulation* sim ){

	emf_report( &sim->emf, BPART, 0 );
	emf_report( &sim->emf, BPART, 1 );
	emf_report( &sim->emf, BPART, 2 );

}
