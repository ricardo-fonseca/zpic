/**
 * ZPIC - em1ds_bnd
 *
 * Laser Wakefield Acceleration
 */

#include "simulation.h"

void sim_init( t_simulation* sim ){

	// Time step
	float dt = 0.019;
	float tmax = 22.8;

	// Simulation box
	int   nx  = 3000;
	float box = 60.0;
	
	// Diagnostic frequency
	int ndump = 50;

    // Initialize particles
	const int n_species = 1;

	// Use 128 particles per cell
	int ppc = 128;

	// Density profile
	t_density density = { .type = SLAB, .start = 20.0, .end = 50.0 };

	t_species* species = (t_species *) malloc( n_species * sizeof( t_species ));
	spec_new( &species[0], "electrons", -1.0, ppc, NULL, NULL, nx, box, dt, &density );


	// Initialize Simulation data
	sim_new( sim, nx, box, dt, tmax, ndump, species, n_species );

	// Add neutralizing background (this must come after sim_new)
	sim_add_neutral_bkg( sim );

	// Add laser pulse (this must come after sim_new)
	t_emf_laser laser = {
		.start = 17.0,
		.fwhm  = 2.0,
		.a0 = 2.0,
		.omega0 = 10.0,
		.polarization = M_PI_2
    };
	sim_add_laser( sim, &laser );

}


void sim_report( t_simulation* sim ){

	// Accelerating field
	emf_report( &sim->emf, EFLD, 0 );

	// Laser field
	emf_report( &sim->emf, EFLD, 2 );
		
	// Charge density - per species
	spec_report( &sim->species[0], CHARGE, NULL, NULL );
    
    // x1u1 phasespace
	const int pha_nx[] = {1024,512};
	const float pha_range[][2] = {{0.0,20.0}, {-2.0,+2.0}};
	spec_report(&sim->species[0], PHASESPACE(X1,U1), pha_nx, pha_range);

	spec_report( &sim->species[0], PARTICLES, NULL, NULL );

}
