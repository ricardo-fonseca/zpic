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

void sim_timings( t_simulation* sim, uint64_t t0, uint64_t t1 ){


	fprintf(stderr, "Time for spec. advance = %f s\n", spec_time());
	fprintf(stderr, "Time for emf   advance = %f s\n", emf_time());
	fprintf(stderr, "Total simulation time  = %f s\n", timer_interval_seconds(t0, t1));
	fprintf(stderr, "\n");


	int npart = 0;
	int i;

	for(i=0; i<sim -> n_species; i++)
		npart += sim -> species[i].np;
	
	float perf = spec_time()/(npart);
	fprintf(stderr, "Particle advance [nsec/part] = %f \n", 1.e9*perf);
	fprintf(stderr, "Particle advance [Mpart/sec] = %f \n", 1.e-6/perf);
}

void sim_new( t_simulation* sim, int nx, float box, float dt, float tmax, int ndump, t_species* species, int n_species ){

	sim -> dt = dt;
	sim -> tmax = tmax;
	sim -> ndump = ndump;

	// Initialize FFT configurations
	fftr_init_cfg( &sim -> fft[FORWARD], nx, FFT_FORWARD, FFTR_CCS);
	fftr_init_cfg( &sim -> fft[BACKWARD], nx, FFT_BACKWARD, FFTR_CCS);
	fftr_init_cfg( &sim -> fft[FORWARD_2], 2*nx, FFT_FORWARD, FFTR_CCS);
	fftr_init_cfg( &sim -> fft[BACKWARD_2], 2*nx, FFT_BACKWARD, FFTR_CCS);


	// Initialize simulation data
	emf_new( &sim -> emf, nx, box, dt, sim -> fft );
	current_new(&sim -> current, nx, box, dt, sim -> fft );
	charge_new(&sim -> charge, nx, box, dt, sim -> fft );

	sim -> n_species = n_species;
	sim -> species = species;

	// Check time step
	float cour = sim->emf.dx;
	if ( dt >= cour ){
		printf("Invalid timestep, courant condition violation, dtmax = %f \n", cour );
		exit(-1);
	}

}

void sim_add_laser( t_simulation* sim,  t_emf_laser* laser ){

	emf_add_laser( &sim->emf, laser );
}

void sim_add_neutral_bkg( t_simulation* sim ){

	unsigned i;
	
	charge_init_neutral_bkg( &sim -> charge );

	for (i = 0; i<sim -> n_species; i++) 
			spec_deposit_charge( &sim -> species[i], sim -> charge.neutral.s );

	charge_update_neutral_bkg( &sim -> charge );

}

void sim_report_energy( t_simulation* sim )
{
	unsigned i;

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

void sim_delete( t_simulation* sim ) {

	unsigned int i;
	for (i = 0; i<sim->n_species; i++) spec_delete( &sim->species[i] );
	
	free( sim->species );

	charge_delete( &sim->charge );
	current_delete( &sim->current );
	emf_delete( &sim->emf );

	fftr_cleanup_cfg( &sim -> fft[FORWARD] );
	fftr_cleanup_cfg( &sim -> fft[BACKWARD] );
	fftr_cleanup_cfg( &sim -> fft[FORWARD_2] );
	fftr_cleanup_cfg( &sim -> fft[BACKWARD_2] );

}

