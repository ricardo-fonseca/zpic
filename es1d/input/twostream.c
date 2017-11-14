/**
 * ZPIC - es1d
 *
 * Two-stream instability
 */

#include <stdlib.h>
#include <math.h>

#include "../simulation.h"

void sim_init( t_simulation* sim ){

	// Time step
	float dt = 0.25;
	float tmax = 100.0;

	// Simulation box
	int   nx  = 120;
	float box = 0.5 * M_PI;

	// Diagnostic frequency
	int ndump = 1;
//	int ndump = 0;

    // Initialize particles
	const int n_species = 2;

	// Use 500 particles per cell
	int ppc = 500;

	// Initial fluid and thermal velocities
	float ufl1 =  0.025;
	float ufl2 = -0.025;

	float uth  = 2e-4;


	t_species* species = (t_species *) malloc( n_species * sizeof( t_species ));
	spec_new( &species[0], "beam1", -1.0, ppc, &ufl1, &uth, nx, box, dt, NULL );
	spec_new( &species[1], "beam2", -1.0, ppc, &ufl2, &uth, nx, box, dt, NULL );

	// Initialize Simulation data
	sim_new( sim, nx, box, dt, tmax, ndump, species, n_species );

}


void sim_report( t_simulation* sim ){

	// Electric field components
	field_report( &sim -> field );

	// Global charge
	charge_report( &sim -> charge );

	// Species particles
	spec_report( &sim->species[0], PARTICLES, NULL, NULL );
	spec_report( &sim->species[1], PARTICLES, NULL, NULL );


}
