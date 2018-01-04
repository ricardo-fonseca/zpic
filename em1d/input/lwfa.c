/**
 * ZPIC - em1d
 *
 * Laser Wakefield Acceleration
 */

#include <stdlib.h>
#include <math.h>

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

	// Density profile
	t_density density = { .type = STEP, .start = 20.0 };

	t_species* species = (t_species *) malloc( n_species * sizeof( t_species ));
	spec_new( &species[0], "electrons", -1.0, ppc, NULL, NULL, nx, box, dt, &density );

	// Initialize Simulation data
	sim_new( sim, nx, box, dt, tmax, ndump, species, n_species );

	// Add laser pulse (this must come after sim_new)
	t_emf_laser laser = {
		.start = 17.0,
		.fwhm  = 2.0,
		.a0 = 2.0,
		.omega0 = 10.0,
		.polarization = M_PI_2
    };
	sim_add_laser( sim, &laser );

	// Set moving window (this must come after sim_new)
	sim_set_moving_window( sim );

	// Set current smoothing (this must come after sim_new)
	t_smooth smooth = {
		.xtype = COMPENSATED,
		.xlevel = 4
	};

	sim_set_smooth( sim, &smooth );

}


void sim_report( t_simulation* sim ){

	// All electric field components
	emf_report( &sim->emf, EFLD, 0 );
	emf_report( &sim->emf, EFLD, 1 );
	emf_report( &sim->emf, EFLD, 2 );

	// Charge density
	spec_report( &sim->species[0], CHARGE, NULL, NULL );

    // x1u1 phasespace
	const int pha_nx[] = {1024,512};
	const float pha_range[][2] = {{0.0,20.0}, {-2.0,+2.0}};
	spec_report(&sim->species[0], PHASESPACE(X1,U1), pha_nx, pha_range);

}
