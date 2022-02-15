#ifndef __SIMULATION__
#define __SIMULATION__

#include <stdint.h>
#include "particles.h"
#include "emf.h"
#include "current.h"

/**
 * @brief EM2D PIC Simulation
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
 * @param sim 	EM1DS simulation
 */
void sim_report( t_simulation* sim );

void sim_iter( t_simulation* sim );
void sim_report_energy( t_simulation* sim );

void sim_new( t_simulation* sim, int nx[], float box[], float dt, float tmax, int ndump, t_species* species, int n_species );
int report( int n, int ndump );
void sim_timings( t_simulation* sim, uint64_t t0, uint64_t t1 );
void sim_add_laser( t_simulation* sim,  t_emf_laser* laser );
void sim_delete( t_simulation* sim );

void sim_set_moving_window( t_simulation* sim );
void sim_set_smooth( t_simulation* sim,  t_smooth* smooth );
void sim_set_ext_fld( t_simulation* sim, t_emf_ext_fld* ext_fld );

#endif
