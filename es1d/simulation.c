#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "simulation.h"
#include "timer.h"


int report( int n, int ndump )
{
	if (ndump > 0) {
		return ! (n % ndump);
	}	else {
		return 0;
	}
}

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


void sim_timings( t_simulation* sim, uint64_t t0, uint64_t t1 ){


	fprintf(stderr, "Time for spec. advance = %f s\n", spec_time());
	fprintf(stderr, "Time for field advance = %f s\n", field_time());
	fprintf(stderr, "Total simulation time  = %f s\n", timer_interval_seconds(t0, t1));
	fprintf(stderr, "\n");


	int npart = 0;

	for(int i=0; i<sim -> n_species; i++)
		npart += sim -> species[i].np;

	if ( npart > 0 ) {
		float perf = spec_time()/(npart);
		fprintf(stderr, "Particle advance [nsec/part] = %f \n", 1.e9*perf);
		fprintf(stderr, "Particle advance [Mpart/sec] = %f \n", 1.e-6/perf);
	}
}

void sim_new( t_simulation* sim, int nx, float box, float dt, float tmax, int ndump, t_species* species, int n_species ){

	sim -> dt = dt;
	sim -> tmax = tmax;
	sim -> ndump = ndump;

	// Initialize FFT configurations
	fftr_init_cfg( &sim -> fft_forward, nx, FFT_FORWARD, FFTR_CCS);
	fftr_init_cfg( &sim -> fft_backward, nx, FFT_BACKWARD, FFTR_CCS);

	// Initialize simulation data
	field_new( &sim -> field, nx, box, dt, &sim->fft_backward );
	charge_new(&sim -> charge, nx, box, dt, &sim->fft_forward );

	sim -> n_species = n_species;
	sim -> species = species;

	// Check time step
	// ... Not sure what kind of test can be done here
}

void sim_add_neutral_bkg( t_simulation* sim ){

	charge_init_neutral_bkg( &sim -> charge );

	for (int i = 0; i<sim -> n_species; i++)
			spec_deposit_charge( &sim -> species[i], sim -> charge.neutral.s );

	charge_update_neutral_bkg( &sim -> charge );

}


void sim_delete( t_simulation* sim ) {

	for (int i = 0; i<sim->n_species; i++) spec_delete( &sim->species[i] );
	
	free( sim->species );

	charge_delete( &sim->charge );
	field_delete( &sim->field );

	fftr_cleanup_cfg( &sim -> fft_backward );
	fftr_cleanup_cfg( &sim -> fft_forward );


}

