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

	int npart = 0;
	int i;

	for(i=0; i<sim -> n_species; i++)
		npart += sim -> species[i].np;

	fprintf(stderr, "Time for spec. advance = %f s\n", spec_time());
	fprintf(stderr, "Time for emf   advance = %f s\n", emf_time());
	fprintf(stderr, "Total simulation time  = %f s\n", timer_interval_seconds(t0, t1));
	fprintf(stderr, "\n");
	
	float perf = spec_time()/(npart);
	fprintf(stderr, "Particle advance [nsec/part] = %f \n", 1.e9*perf);
	fprintf(stderr, "Particle advance [Mpart/sec] = %f \n", 1.e-6/perf);
}

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

void sim_add_laser( t_simulation* sim,  t_emf_laser* laser ){

	emf_add_laser( &sim->emf, laser );
}

void sim_set_smooth( t_simulation* sim,  t_smooth* smooth ){
    
    if ( (smooth -> xtype != none) && (smooth -> xlevel <= 0) ) {
    	printf("Invalid smooth level along x direction\n");
    	exit(-1);
    }
	sim -> current.smooth = *smooth;
}

void sim_set_moving_window( t_simulation* sim ){

	// Set moving window flag and disable boundary conditions
	// for EM fields
	sim -> emf.moving_window = 1;
    sim -> emf.bc_type = EMF_BC_NONE;

	// Disable boundary conditions for electric current
	// No specific code for 
	sim -> current.bc_type = CURRENT_BC_NONE;

	// Set moving window flag for all species
	for(int i=0; i<sim -> n_species; i++)
		sim -> species[i].moving_window = 1;

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

	int i;
	for (i = 0; i<sim->n_species; i++) spec_delete( &sim->species[i] );
	
	free( sim->species );

	current_delete( &sim->current );
	emf_delete( &sim->emf );

}

