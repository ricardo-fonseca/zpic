/**
 * @file simulation.c
 * @author Ricardo Fonseca
 * @brief EM2DS Simulation
 * @version 0.2
 * @date 2022-02-18
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
 * @brief Checks if the `sim_report()` function should be called 
 * 
 * @param n 		Current iteration
 * @param ndump 	Diagnostic frequency (number of iterations between diagnostic dumps)
 * @return 			1 if `sim_report()` function should be called, 0 otherwise
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
 * 1. Zeroing current/charge density
 * 2. Advancing particle species and depositing electric current/charge
 * 3. Updating electric current/charge boundary
 * 4. Advancing the EM fields
 * 
 * @param sim 	EM2DS Simulation
 */
void sim_iter( t_simulation* sim ) {
	// Advance particles and deposit current
	current_zero( &sim -> current );
	charge_zero( &sim -> charge );
	for (int i = 0; i<sim -> n_species; i++)
		spec_advance(&sim -> species[i], &sim -> emf, &sim -> charge, &sim -> current );

	// Update charge and current boundary conditions and get fourier transforms
	current_update( &sim -> current );
	charge_update( &sim -> charge );

	// Advance EM fields
	emf_advance( &sim -> emf, &sim -> charge, &sim -> current );

}

/**
 * @brief Prints out report on simulation timings
 * 
 * @param sim 	EM2DS Simulaiton
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
 * @param sim 			EM2DS Simulation
 * @param nx 			Number of grid points [x,y]
 * @param box 			Simulation box size in simulation units [x,y]
 * @param dt 			Simulation time step in simulation units
 * @param tmax 			Final simulation time
 * @param ndump 		Diagnostic frequency (`report` function will be called every ndump iterations)
 * 						set to 0 to disable diagnostic reports
 * @param species 		Array of particle species, may be NULL (no particles)
 * @param n_species 	Number of particle specis
 */
void sim_new( t_simulation* sim, const int nx[], float box[], float dt, float tmax, int ndump, t_species* species, int n_species ){

	sim -> dt = dt;
	sim -> tmax = tmax;
	sim -> ndump = ndump;

	emf_new( &sim -> emf, nx, box, dt );
	current_new(&sim -> current, nx, box, dt);
	charge_new(&sim -> charge, nx, box, dt );

	sim -> n_species = n_species;
	sim -> species = species;

	// Check time step
	float cour = sqrtf( 1.0f/( 1.0f/(sim->emf.dx[0]*sim->emf.dx[0]) + 1.0f/(sim->emf.dx[1]*sim->emf.dx[1]) ) );
	if ( dt >= cour ){
		printf("Invalid timestep, courant condition violation, dtmax = %f \n", cour );
		exit(-1);
	}

}

/**
 * @brief Adds laser pulse to simulation
 * 
 * This routine is just used to avoid users accessing the simulation.emf object directly,
 * see `emf_add_laser()` for details.
 * 
 * @param sim 		EM2DS simulation
 * @param laser 	Laser parameters
 */
void sim_add_laser( t_simulation* sim,  t_emf_laser* laser ){

	emf_add_laser( &sim->emf, laser );
}

/**
 * @brief Adds neutralizing background to the simulation
 * 
 * The neutralizing background is calculated by depositing the charge from
 * all particle species
 * 
 * @param sim 	EM2DS Simulation
 */
void sim_add_neutral_bkg( t_simulation* sim ){

	charge_init_neutral_bkg( &sim -> charge );

	for (int i = 0; i<sim -> n_species; i++) 
			spec_deposit_charge( &sim -> species[i], sim -> charge.neutral.s );

	charge_update_neutral_bkg( &sim -> charge );

}

/**
 * @brief Print report on simulation energy (fields/particles/total)
 * 
 * @param sim 	EM2DS Simulation
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
 * @param sim 		EM2DS Simulation
 */
void sim_delete( t_simulation* sim ) {

	int i;
	for (i = 0; i<sim->n_species; i++) spec_delete( &sim->species[i] );

	free( sim->species );

	charge_delete( &sim->charge );
	current_delete( &sim->current );
	emf_delete( &sim->emf );

}

