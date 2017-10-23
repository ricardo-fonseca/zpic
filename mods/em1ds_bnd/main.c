/*
Copyright (C) 2017 Instituto Superior Tecnico

This file is part of the ZPIC Educational code suite

The ZPIC Educational code suite is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

The ZPIC Educational code suite is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with the ZPIC Educational code suite. If not, see <http://www.gnu.org/licenses/>.
*/

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
#include "input/lwfa.c"

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
			spec_advance( &sim.species[i], &sim.emf, &sim.charge, &sim.current );
		
		// Update cahrge and current boundary conditions and get fourier transforms
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
