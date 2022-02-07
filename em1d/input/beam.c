/**
 * ZPIC - em1d
 *
 * Initializing fluid velocity
 */

#include <stdlib.h>
#include "../simulation.h"

void sim_init( t_simulation* sim ){

	// Time step
	float dt = 0.019;
	float tmax = 22.8;

	// Simulation box
	int   nx  = 1000;
	float box = 20.0;

	// Diagnostic frequency
	int ndump = 50;

    // Initialize particles
	const int n_species = 1;

	// Use 128 particles per cell
	int ppc = 128;

	// Initial fluid and thermal velocities
	float ufl[] = { 0.2, 0.0, 0.0 };
	float uth[] = { 0.0, 0.0, 0.0 };



	t_species* species = (t_species *) malloc( n_species * sizeof( t_species ));
	spec_new( &species[0], "electrons", -1.0, ppc, ufl, uth, nx, box, dt, NULL );

	// Initialize Simulation data
	sim_new( sim, nx, box, dt, tmax, ndump, species, n_species );

}

void sim_report( t_simulation* sim ){

	// All electric field components
	emf_report( &sim->emf, EFLD, 0 );
	emf_report( &sim->emf, EFLD, 1 );
	emf_report( &sim->emf, EFLD, 2 );

	// Charge density
	spec_report( &sim->species[0], CHARGE, NULL, NULL );

    // x1u1 phasespace
	const int pha_nx[] = {1000,512};
	const float pha_range[][2] = {{0.0,20.0}, {0.0,1.5}};
	spec_report(&sim->species[0], PHASESPACE(X1,U1), pha_nx, pha_range);

}
