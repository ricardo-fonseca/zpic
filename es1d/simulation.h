#ifndef __SIMULATION__
#define __SIMULATION__

#include <stdint.h>
#include <stdlib.h>
#include <math.h>

#include "particles.h"
#include "field.h"
#include "charge.h"

/**
 * @brief ES1D PIC 
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
	t_field field;			///< Electric field
	t_charge charge;		///< Electric charge density

	// FFT
	t_fftr_cfg fft_forward;		///< Forward (real to complex) FFT configuration
	t_fftr_cfg fft_backward;	///< Backward (complex to real) FFT configuration

} t_simulation;


void sim_init( t_simulation* sim );
void sim_report( t_simulation* sim );

void sim_iter( t_simulation* sim );
void sim_report_energy( t_simulation* sim );

void sim_add_neutral_bkg( t_simulation* sim );

void sim_new( t_simulation* sim, int nx, float box, float dt, float tmax, int ndump, t_species* species, int n_species );
int report( int n, int ndump );
void sim_timings( t_simulation* sim, uint64_t t0, uint64_t t1 );
void sim_delete( t_simulation* sim );

#endif
