/**
 * @file simulation.c
 * @author Ricardo Fonseca
 * @brief ES1D Simulation
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
 * 1. Zeroing charge density
 * 2. Advancing particle species and depositing electric charge
 * 3. Updating electric charge boundary
 * 4. Advancing the electric field
 * 
 * @param sim 	ES1D Simulation
 */
void sim_iter( t_simulation* sim ) {
	// Advance particles and deposit charge
	charge_zero( &sim -> charge );
	for (int i = 0; i < sim -> n_species; i++)
		spec_advance(&sim -> species[i], &sim -> field, &sim  -> charge );

	// Update charge boundary conditions and get fourier transforms
	charge_update( &sim -> charge );

	// Advance E field
	field_advance( &sim -> field, &sim -> charge );
}

/**
 * @brief Prints out report on simulation timings
 * 
 * @param sim 	ES1D Simulaiton
 * @param t0 	Simulation start time (ticks)
 * @param t1 	Simulation end time (ticks)
 */
void sim_timings( t_simulation* sim, uint64_t t0, uint64_t t1 ){

	fprintf(stderr, "Time for spec. advance = %f s\n", spec_time());
	fprintf(stderr, "Time for field advance = %f s\n", field_time());
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
 * @param sim 			ES1D Simulation
 * @param nx 			Number of grid points
 * @param box 			Simulation box size in simulation units
 * @param dt 			Simulation time step in simulation units
 * @param tmax 			Final simulation time
 * @param ndump 		Diagnostic frequency (`report` function will be called every ndump iterations)
 * 						set to 0 to disable diagnostic reports
 * @param species 		Array of particle species, may be NULL (no particles)
 * @param n_species 	Number of particle specis
 */
void sim_new( t_simulation* sim, int nx, float box, float dt, float tmax, int ndump, t_species* species, int n_species ){

	sim -> dt = dt;
	sim -> tmax = tmax;
	sim -> ndump = ndump;

	// Initialize FFT configurations
	fftr_init_cfg( &sim -> fft_forward, nx, FFT_FORWARD );
	fftr_init_cfg( &sim -> fft_backward, nx, FFT_BACKWARD );

	// Initialize simulation data
	field_new( &sim -> field, nx, box, dt, &sim->fft_backward );
	charge_new(&sim -> charge, nx, box, dt, &sim->fft_forward );

	sim -> n_species = n_species;
	sim -> species = species;

	// Check time step
	// ... Not sure what kind of test can be done here
}

/**
 * @brief Adds neutralizing background to the simulation
 * 
 * The neutralizing background is calculated by depositing the charge from
 * all particle species
 * 
 * @param sim 	ES1D Simulation
 */
void sim_add_neutral_bkg( t_simulation* sim ){

	charge_init_neutral_bkg( &sim -> charge );

	for (int i = 0; i<sim -> n_species; i++)
			spec_deposit_charge( &sim -> species[i], sim -> charge.neutral.s );

	charge_update_neutral_bkg( &sim -> charge );

}

/**
 * @brief Frees dynamic memory 
 * 
 * @param sim 		ES1D Simulation
 */
void sim_delete( t_simulation* sim ) {

	for (int i = 0; i<sim->n_species; i++) spec_delete( &sim->species[i] );
	
	free( sim->species );

	charge_delete( &sim->charge );
	field_delete( &sim->field );

	fftr_cleanup_cfg( &sim -> fft_backward );
	fftr_cleanup_cfg( &sim -> fft_forward );


}

