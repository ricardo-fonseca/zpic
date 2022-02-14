#ifndef __SIMULATION__
#define __SIMULATION__

#include <stdint.h>
#include <stdlib.h>
#include <math.h>

#include "filter.h"
#include "particles.h"
#include "emf.h"
#include "current.h"

/**
 * @brief EM1DS PIC Simulation
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
	t_charge charge;		///< Electric charge density

	// FFT
	t_fftr_cfg fft_forward;		///< Forward (real to complex) FFT configuration
	t_fftr_cfg fft_backward;	///< Backward (complex to real) FFT configuration

	/// Spectral filtering
	t_filter filter;

} t_simulation;


/**
 * @brief Initializes simulation parameters
 * 
 * This routine __MUST__ be supplied by the user, see the `input` directory
 * for examples
 * 
 * @param sim	EM1DS simulation 
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
 * @param sim 	EM1DS Simulation
 */
void sim_report( t_simulation* sim );

/**
 * @brief Advance simulation 1 iteration
 * 
 * A complete iteration consists of:
 * 1. Zeroing current/charge density
 * 2. Advancing particle species and depositing electric current/charge
 * 3. Updating electric current/charge boundary
 * 4. Advancing the EM fields
 * 
 * @param sim 	EM1DS Simulation
 */
void sim_iter( t_simulation* sim );

/**
 * @brief Print report on simulation energy (fields/particles/total)
 * 
 * @param sim 	EM1DS Simulation
 */
void sim_report_energy( t_simulation* sim );

/**
 * @brief Adds laser pulse to simulation
 * 
 * This routine is just used to avoid users accessing the simulation.emf object directly,
 * see `emf_add_laser()` for details.
 * 
 * @param sim 		EM1DS simulation
 * @param laser 	Laser parameters
 */
void sim_add_laser( t_simulation* sim,  t_emf_laser* laser );

/**
 * @brief Adds neutralizing background to the simulation
 * 
 * The neutralizing background is calculated by depositing the charge from
 * all particle species
 * 
 * @param sim 	EM1DS Simulation
 */
void sim_add_neutral_bkg( t_simulation* sim );

/**
 * @brief Sets filtering parameters for the simulation
 * 
 * @param sim 		EM1DS Simulartion
 * @param type 		Filter type {FILTER_NONE, FILTER_GAUSS, FILTER_SHARP}
 * @param ck 		Cutoff wavenumber for the filter (filter response = 1/2)
 * @return 			0 on success, -1 on error
 */
int sim_filter_set( t_simulation* const sim, enum filter_type const type, float const ck );

/**
 * @brief Initialize simulation object
 * 
 * @param sim 			EM1DS Simulation
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
 * @param sim 	EM1DS Simulaiton
 * @param t0 	Simulation start time (ticks)
 * @param t1 	Simulation end time (ticks)
 */
void sim_timings( t_simulation* sim, uint64_t t0, uint64_t t1 );

/**
 * @brief Frees dynamic memory 
 * 
 * @param sim 		EM1DS Simulation
 */
void sim_delete( t_simulation* sim );

/**
 * @brief Checks if the `sim_report()` function should be called 
 * 
 * @param n 		Current iteration
 * @param ndump 	Diagnostic frequency (number of iterations between diagnostic dumps)
 * @return 			1 if `sim_report()` function should be called, 0 otherwise
 */
int report( int n, int ndump );

#endif
