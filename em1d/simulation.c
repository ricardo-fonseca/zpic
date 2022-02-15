/**
 * @file simulation.c
 * @author Ricardo Fonseca
 * @brief EM1D Simulation
 * @version 0.2
 * @date 2022-02-04
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "simulation.h"
#include "timer.h"

/**
 * @brief Checks if there should be a report at this timestep
 * 
 * @param n 		Current iteration
 * @param ndump 	Diagnostic frequency (number of iterations between diagnostic dumps)
 * @return 			1 if it is time to write a report, 0 otherwise
 */
int report( int n, int ndump )
{
	if (ndump > 0) {
		return ! (n % ndump);
	}	else {
		return 0;
	}
}

/**
 * @brief Advance simulation 1 iteration
 * 
 * A complete iteration consists of:
 * 1. Zeroing current density
 * 2. Advancing particle species and depositing electric current
 * 3. Updating electric current boundary
 * 4. Advancing the EM fields
 * 
 * @param sim 	EM1D Simulation
 */
void sim_iter( t_simulation* sim ) {
	// Advance particles and deposit current
	current_zero( &sim -> current );
	for (int i = 0; i<sim -> n_species; i++)
		spec_advance(&sim -> species[i], &sim -> emf, &sim -> current );

	// Update current boundary conditions and advance iteration
	current_update( &sim -> current );

	// Advance EM fields
	emf_advance( &sim -> emf, &sim -> current );
}

/**
 * @brief Prints out report on simulation timings
 * 
 * @param sim 	EM1D Simulaiton
 * @param t0 	Simulation start time (ticks)
 * @param t1 	Simulation end time (ticks)
 */
void sim_timings( t_simulation* sim, uint64_t t0, uint64_t t1 ){

	fprintf(stderr, "Time for spec. advance = %f s\n", spec_time());
	fprintf(stderr, "Time for emf   advance = %f s\n", emf_time());
	fprintf(stderr, "Total simulation time  = %f s\n", timer_interval_seconds(t0, t1));
	fprintf(stderr, "\n");

	double perf = spec_perf();
	if ( perf > 0 ) {
		fprintf(stderr, "Particle advance [nsec/part] = %f \n", 1.e9*perf);
		fprintf(stderr, "Particle advance [Mpart/sec] = %f \n", 1.e-6/perf);
	}
}

/**
 * @brief Initialize simulation object
 * 
 * @param sim 			EM1D Simulation
 * @param nx 			Number of grid points
 * @param box 			Simulation box size in simulation units
 * @param dt 			Simulation time step in simulation units
 * @param tmax 			Final simulation time
 * @param ndump 		Diagnostic frequency (`sim_report()` function will be called every ndump iterations)
 * 						set to 0 to disable diagnostic reports
 * @param species 		Array of particle species, may be NULL (no particles)
 * @param n_species 	Number of particle specis
 */
void sim_new( t_simulation* sim, int nx, float box, float dt, float tmax, int ndump, t_species* species, int n_species ){

	sim -> dt = dt;
	sim -> tmax = tmax;
	sim -> ndump = ndump;

	emf_new( &sim -> emf, nx, box, dt );
	current_new(&sim -> current, nx, box, dt);

	sim -> n_species = n_species;
	sim -> species = species;

	// Check time step
	float cour = sim->emf.dx;
	if ( dt >= cour ){
		fprintf(stderr, "Invalid timestep, courant condition violation, dtmax = %f \n", cour );
		exit(-1);
	}

}

/**
 * @brief Adds laser pulse to simulation
 * 
 * This routine is just used to avoid users accessing the simulation.emf object directly,
 * see `emf_add_laser()` for details.
 * 
 * @param sim 		EM1D simulation
 * @param laser 	Laser parameters
 */
void sim_add_laser( t_simulation* sim,  t_emf_laser* laser ){
	emf_add_laser( &sim->emf, laser );
}

/**
 * @brief Sets external EM fields for the simulation
 * 
 * This routine is just used to avoid users accessing the simulation.emf object directly,
 * see `emf_set_ext_fld()` for details.
 * 
 * @param sim 		EM1D simulation
 * @param ext_fld 	External EM fields parameters
 */
void sim_set_ext_fld( t_simulation* sim, t_emf_ext_fld* ext_fld ){
	emf_set_ext_fld( &sim->emf, ext_fld );
}

/**
 * @brief Sets electric current digital filtering (smoothing) parameters
 * 
 * @param sim 		EM1D Simulation
 * @param smooth 	Digital filtering parameters
 */
void sim_set_smooth( t_simulation* sim,  t_smooth* smooth ){

    if ( (smooth -> xtype != NONE) && (smooth -> xlevel <= 0) ) {
    	printf("Invalid smooth level along x direction\n");
    	exit(-1);
    }
	sim -> current.smooth = *smooth;
}

/**
 * @brief Sets a moving window algorithm for the simulation
 * 
 * @param sim 
 */
void sim_set_moving_window( t_simulation* sim ){

	// Set moving window flag and disable boundary conditions
	// for EM fields
	sim -> emf.moving_window = 1;
    sim -> emf.bc_type = EMF_BC_NONE;

	// Disable boundary conditions for electric current
	sim -> current.bc_type = CURRENT_BC_NONE;

	// Set moving window flag for all species
	for(int i=0; i<sim -> n_species; i++)
		sim -> species[i].moving_window = 1;
}

/**
 * @brief Print report on simulation energy (fields/particles/total)
 * 
 * @param sim 	EM1D Simulation
 */
void sim_report_energy( t_simulation* sim )
{
	int i;

	double emf_energy[6];
	double part_energy[ sim -> n_species ];

	emf_get_energy( &sim -> emf, emf_energy );
	double tot_emf = emf_energy[0];
	for( i = 0; i < 6; i++ ){
		tot_emf += emf_energy[i];
	}

	double tot_part = 0;
	for( i = 0; i < sim -> n_species; i++ ){
		part_energy[i] = sim -> species[i].energy;
		tot_part += part_energy[i];
	}

	printf("Energy (fields | particles | total) = %e %e %e\n",
		tot_emf, tot_part, tot_emf+tot_part);

}

/**
 * @brief Frees dynamic memory 
 * 
 * @param sim 		EM1D Simulation
 */
void sim_delete( t_simulation* sim ) {

	for (int i = 0; i<sim->n_species; i++) spec_delete( &sim->species[i] );

	free( sim->species );

	current_delete( &sim->current );
	emf_delete( &sim->emf );

}
