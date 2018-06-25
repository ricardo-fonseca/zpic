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

void sim_timings( t_simulation* sim, uint64_t t0, uint64_t t1 ){

	fprintf(stderr, "Time for spec. advance = %f s\n", spec_time());
	fprintf(stderr, "Time for emf   advance = %f s\n", emf_time());
	fprintf(stderr, "Total simulation time  = %f s\n", timer_interval_seconds(t0, t1));
	fprintf(stderr, "\n");

	if (spec_time()>0) {
		double perf = spec_perf();
		fprintf(stderr, "Particle advance [nsec/part] = %f \n", 1.e9*perf);
		fprintf(stderr, "Particle advance [Mpart/sec] = %f \n", 1.e-6/perf);
	}
}

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

void sim_add_laser( t_simulation* sim,  t_emf_laser* laser ){

	emf_add_laser( &sim->emf, laser );
}

void sim_add_neutral_bkg( t_simulation* sim ){

	charge_init_neutral_bkg( &sim -> charge );

	for (int i = 0; i<sim -> n_species; i++) 
			spec_deposit_charge( &sim -> species[i], sim -> charge.neutral.s );

	charge_update_neutral_bkg( &sim -> charge );

}


void sim_delete( t_simulation* sim ) {

	int i;
	for (i = 0; i<sim->n_species; i++) spec_delete( &sim->species[i] );

	free( sim->species );

	charge_delete( &sim->charge );
	current_delete( &sim->current );
	emf_delete( &sim->emf );

}

