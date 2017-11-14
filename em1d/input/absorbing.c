/**
 * ZPIC - em1d
 *
 * Demonstration of the absorbing boundary conditions
 */

#include <stdlib.h>
#include <math.h>

#include "../simulation.h"

void sim_init( t_simulation* sim ){

	// Time step
	float dt = 0.019;
//	float tmax = 22.8;
	float tmax = 40.0;

	// Simulation box
	int   nx  = 1000;
	float box = 20.0;

	// Diagnostic frequency
	int ndump = 50;

    // Initialize particles
	const int n_species = 0;

	// Initialize Simulation data
	sim_new( sim, nx, box, dt, tmax, ndump, NULL, n_species );

	// Add laser pulse (this must come after sim_new)
	t_emf_laser laser = {
		.start = 17.0,
		.fwhm  = 2.0,
		.a0 = 2.0,
		.omega0 = 10.0,
		.polarization = M_PI_2
    };
	sim_add_laser( sim, &laser );

	sim -> emf.bc_type = EMF_BC_OPEN;

}


void sim_report( t_simulation* sim ){

	// All electric field components
	emf_report( &sim->emf, EFLD, 0 );
	emf_report( &sim->emf, EFLD, 1 );
	emf_report( &sim->emf, EFLD, 2 );


	emf_report( &sim->emf, BFLD, 0 );
	emf_report( &sim->emf, BFLD, 1 );
	emf_report( &sim->emf, BFLD, 2 );

}
