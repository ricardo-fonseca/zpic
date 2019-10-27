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
#include <unistd.h>
#ifdef __linux__
#   include <linux/limits.h>
#else
#   include <limits.h>
#endif

#ifdef _MPI_
#include <mpi.h>
#endif

#include "zpic.h"
#include "simulation.h"
#include "emf.h"
#include "current.h"
#include "particles.h"
#include "timer.h"
#include "zdf.h"

// Input file
#include "input/lwfa.c"


/**
 * Run a simulation for the selected parameter set
 * @param opt_idx Parameter set index
 */
void run_sim( int opt_idx ) {

	// Initialize simulation
	t_simulation sim;
	sim_init( &sim, opt_idx );
	
	// Get working directory
	char root[PATH_MAX];
	getcwd(root, PATH_MAX);

	// Create output directory (if needed) and move there
	char outpath[NAME_MAX];
	snprintf(outpath, NAME_MAX, "output.%d", opt_idx);
	create_path( outpath );
	chdir( outpath );

    // Run simulation
	int n, i;
	float t;
	
	printf("Starting simulation %d...\n", opt_idx);
		
	for (n=0,t=0.0; t<=sim.tmax; n++, t=n*sim.dt) {
		
		// Don't print progress

		if ( report ( n , sim.ndump ) ) 
			sim_report( &sim );
		
		// Advance particles and deposit current
		current_zero( &sim.current );
		for (i = 0; i<sim.n_species; i++) 
			spec_advance(&sim.species[i], &sim.emf, &sim.current );
		
		// Update current boundary conditions and advance iteration
		current_update( &sim.current );
		
		// Advance EM fields
		emf_advance( &sim.emf, &sim.current );		
	}

	printf("Simulation %d ended.\n", opt_idx);

    // Return to original working directory
    chdir(root);
	
    // Cleanup data
    sim_delete( &sim );
}


int main (int argc, char * argv[]) {

#ifdef _MPI_

	int myrank, size;

	// Initialize MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	if ( size >= n_scan ) {

		// If more processes than options, run 1 option per process and keep
		// excess processes idle
		if ( myrank < n_scan ) run_sim(myrank);
		else printf("Process %d is idle.\n", myrank);

	} else {

		// Divide parameters as evenly as possible accross processes
		int st = n_scan / size ;
		int rm = n_scan % size ;
		int i0, i1;

		if ( myrank < rm ) {
			i0 = myrank * (st + 1);
			i1 = i0 + st;
		} else {
			i0 = myrank * st + rm;
			i1 = i0 + (st-1); 
		}

		for( int i = i0; i <= i1; i++ ) run_sim(i);
	}

	// Finalize MPI
	MPI_Finalize();

#else
    // Run options sequentially
	for (int i = 0; i < n_scan; i++) {
		run_sim(i);
	}
#endif

	return 0;
}
