/**
 * ZPIC - em2ds
 *
 * Launch a laser pulse
 */

#include <stdlib.h>
#include <math.h>

#include "../simulation.h"

void sim_init( t_simulation* sim ){

	// Time step
	float dt = 0.014;
	float tmax = 20.0;

	// Simulation box
	int   nx[2]  = { 1000, 128 };
	float box[2] = { 20.0, 25.6 };

	// Diagnostic frequency
	int ndump = 100;


	// Initialize Simulation data
	sim_new( sim, nx, box, dt, tmax, ndump, NULL, 0 );

	// Add laser pulse (this must come after sim_new)
	t_emf_laser laser = {
		.type = GAUSSIAN,
		.start = 17.0,
		.fwhm  = 2.0,
		.a0 = 2.0,
		.omega0 = 10.0,
		.W0 = 4.0,
		.focus = 20.0,
		.axis = 12.8,
		.polarization = M_PI_2
    };
	sim_add_laser( sim, &laser );

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
