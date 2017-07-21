#include <stdio.h>
#include <stdlib.h>

#include <math.h>

#include "zpic.h"
#include "simulation.h"
#include "emf.h"
#include "current.h"
#include "charge.h"
#include "particles.h"
#include "timer.h"

// Include Simulation parameters here
//#include "lwfa.c"
//#include "twostream.c"
#include "density.c"

int main (int argc, const char * argv[]) {
    
	// Initialize simulation
	t_simulation sim;
	sim_init( &sim );
	
    // Run simulation
	int n, i;
	float t;
	
	printf("Starting simulation ...\n\n");
	
	uint64_t t0,t1;
	t0 = timer_ticks();
	
	for (n=0,t=0.0; t<=sim.tmax; n++, t=n*sim.dt) {
		printf("n = %i, t = %f\n",n,t);

		if ( report ( n , sim.ndump ) ) 
			sim_report( &sim );
		
		// Advance particles and deposit charge and current
		current_zero( &sim.current );
		charge_zero( &sim.charge );
		for (i = 0; i<sim.n_species; i++) 
			spec_advance(&sim.species[i], &sim.emf, &sim.charge, &sim.current );
		
		// Update charge and current boundary conditions and get fourier transforms
		current_update( &sim.current );
		charge_update( &sim.charge );
		
		// Advance EM fields
		emf_advance( &sim.emf, &sim.charge, &sim.current );		
	}

	t1 = timer_ticks();
	fprintf(stderr, "\nSimulation ended.\n\n");

	// Simulation times
    sim_timings( &sim, t0, t1 );
	
    // Cleanup data
    sim_delete( &sim );
    
	return 0;
}

