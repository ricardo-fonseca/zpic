#include <stdio.h>
#include <stdlib.h>

#include <math.h>

#include "zpic.h"
#include "simulation.h"
#include "emf.h"
#include "current.h"
#include "particles.h"
#include "timer.h"

// Include Simulation parameters here
//#include "weibel.c"
#include "lwfa.c"

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
			sim_report( &sim.emf, &sim.current, sim.species );
		
		// Advance particles and deposit current
		current_zero( &sim.current );
		for (i = 0; i<sim.n_species; i++) 
			spec_advance(&sim.species[i], &sim.emf, &sim.current );
		
		// Update current boundary conditions and advance iteration
		current_update( &sim.current );
		
		// Advance EM fields
		emf_advance( &sim.emf, &sim.current );		
	}

	t1 = timer_ticks();
	fprintf(stderr, "\nSimulation ended.\n\n");

	// Simulation times
    sim_timings( &sim, t0, t1 );
	
    // Cleanup data
    sim_delete( &sim );
    
	return 0;
}
