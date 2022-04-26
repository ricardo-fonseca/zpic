#ifndef __SIMULATION__
#define __SIMULATION__

#include <stdint.h>
#include "particles.h"
#include "emf.h"
#include "current.h"

/**
 * @brief EM1D PIC Simulation
 * 
 */
typedef struct Simulation {

	// Time
	float dt;		///< Time step
	float tmax;		///< Final simulation time

	int ndump;		///< Diagnostic frequency

	// Simulation data
	int n_species;			///< Number of particle species
	t_species* species;		///< Particle species
	t_emf emf;				///< EM fields
	t_current current;		///< Electric current density

	int moving_window;		///< Use moving window

} t_simulation;


/**
 * @brief Initializes simulation parameters
 * 
 * This routine __MUST__ be supplied by the user, see the `input` directory
 * for examples
 * 
 * @param sim	EM1D simulation 
 */
void sim_init( t_simulation* sim );

/**
 * @brief Saves diagnostic information
 *
 * This routine __MUST__ be supplied by the user, see the `input` directory
 * for examples
 * 
 * This routine will be called every `ndump` iterations
 * 
 * @param sim 	EM1D simulation
 */
void sim_report( t_simulation* sim );

/**
 * @brief Advance simulation 1 iteration
 * 
 * @param sim 	EM1D Simulation
 */
void sim_iter( t_simulation* sim );

/**
 * @brief Print report on simulation energy (fields/particles/total)
 * 
 * @param sim 	EM1D Simulation
 */
void sim_report_energy( t_simulation* sim );

/**
 * @brief Initialize simulation object
 * 
 * @param sim 			EM1D Simulation
 * @param nx 			Number of grid points
 * @param box 			Simulation box size in simulation units
 * @param dt 			Simulation time step in simulation units
 * @param tmax 			Final simulation time
 * @param ndump 		Diagnostic frequency (`report` function will be called every ndump iterations)
 * 						set to 0 to disable diagnostic reports
 * @param species 		Array of particle species, may be NULL (no particles)
 * @param n_species 	Number of particle specis
 */
void sim_new( t_simulation* sim, int nx, float box, float dt, float tmax, int ndump, t_species* species, int n_species );

/**
 * @brief Prints out report on simulation timings
 * 
 * @param sim 	EM1D Simulaiton
 * @param t0 	Simulation start time (ticks)
 * @param t1 	Simulation end time (ticks)
 */
void sim_timings( t_simulation* sim, uint64_t t0, uint64_t t1 );

/**
 * @brief Adds laser pulse to simulation
 * 
 * @param sim 		EM1D simulation
 * @param laser 	Laser parameters
 */
void sim_add_laser( t_simulation* sim,  t_emf_laser* laser );

/**
 * @brief Frees dynamic memory 
 * 
 * @param sim 		EM1D Simulation
 */
void sim_delete( t_simulation* sim );

/**
 * @brief Sets a moving window algorithm for the simulation
 * 
 * @param sim 
 */
void sim_set_moving_window( t_simulation* sim );

/**
 * @brief Sets electric current digital filtering (smoothing) parameters
 * 
 * @param sim 		EM1D Simulation
 * @param smooth 	Digital filtering parameters
 */
void sim_set_smooth( t_simulation* sim,  t_smooth* smooth );

/**
 * @brief Sets external EM fields for the simulation
 * 
 * @param sim 		EM1D simulation
 * @param ext_fld 	External EM fields parameters
 */
void sim_set_ext_fld( t_simulation* sim, t_emf_ext_fld* ext_fld );

/**
 * @brief Checks if the `sim_report()` function should be called 
 * 
 * @param n 		Current iteration
 * @param ndump 	Diagnostic frequency (number of iterations between diagnostic dumps)
 * @return 			1 if `sim_report()` function should be called, 0 otherwise
 */
int report( int n, int ndump );

#endif
