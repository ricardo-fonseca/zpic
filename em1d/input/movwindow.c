/**
 * ZPIC - em1d
 *
 * Moving simulation window with different density profiles
 */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "../simulation.h"

void sim_init( t_simulation* sim ){
	// Time step
	float dt = 0.08;
	float tmax = 300.00;

	// Simulation box
	int   nx  = 1024;
	float box = 82.0;

	// Diagnostic frequency
	int ndump = 50;

    // Initialize particles
	const int n_species = 1;

	// Use 128 particles per cell
	int ppc = 256;

	// Density profile
//	t_density density = { .type = UNIFORM };
	t_density density = { .type = STEP, .start = 80 };
//	t_density density = { .type = SLAB, .start = 17.5, .end = 22.5 };
//	t_density density = { .type = RAMP, .start = 275, .end = 400, .ramp = {0.0, 150.} };

	t_species* species = (t_species *) malloc( n_species * sizeof( t_species ));
	spec_new( &species[0], "electrons", -1.0, ppc, NULL, NULL, nx, box, dt, &density );

	// Initialize Simulation data
	sim_new( sim, nx, box, dt, tmax, ndump, species, n_species );

	// Add laser pulse (this must come after sim_new)
	t_emf_laser laser = {
	  .start =50.0,
	  .fwhm  = 15,
	  .a0 = 0.1,
	  .omega0 = 10,
	  .polarization = M_PI_2
	};
	sim_add_laser( sim, &laser );

    t_emf_laser laser2 = {
        .start =50.0,
        .fwhm  = 15,
        .a0 = 0.01,
        .omega0 = 11,
        .polarization = M_PI_2
    };
    sim_add_laser( sim, &laser2 );


    t_emf_laser laser3 = {
        .start =50.0,
        .fwhm  = 15,
        .a0 = 0.01,
        .omega0 = 9,
        .polarization = M_PI_2
    };
    sim_add_laser( sim, &laser3 );


	// Set moving window (this must come after sim_new)
	sim_set_moving_window( sim );

}


void sim_report( t_simulation* sim ){

  	// All electric field components
  //	emf_report( &sim->emf, EFLD, 0 );
	emf_report( &sim->emf, EFLD, 2 );
	//	emf_report( &sim->emf, EFLD, 2 );

	// Charge density
	spec_report( &sim->species[0], CHARGE, NULL, NULL );

	// RAW dump
	//	spec_report( &sim->species[0], PARTICLES, NULL, NULL );

}
