/*
 *  particles.h
 *  zpic
 *
 *  Created by Ricardo Fonseca on 11/8/10.
 *  Copyright 2010 Centro de Física dos Plasmas. All rights reserved.
 *
 */

#ifndef __PARTICLES__
#define __PARTICLES__

#include "zpic.h"
#include "emf.h"
#include "current.h"

#define MAX_SPNAME_LEN 32

typedef struct {
    int ix, iy;
    t_part_data x, y;
    t_part_data ux, uy, uz;
} t_part;

enum density_type {UNIFORM, EMPTY, STEP, SLAB, CUSTOM};

typedef struct {

    float n;				// reference density (defaults to 1.0, multiplies density profile)

    enum density_type type;		// Density profile type
    float start, end;		// Position of the plasma start/end, in simulation units

    // Custom density profile parameters
    
    // Pointer to custom density function along x
    float (*custom_x)(float, void*);
    // Pointer to additional data to be passed to the custom_x function
    void *custom_data_x;
    // Pointer to custom density function along y
    float (*custom_y)(float, void*);
    // Pointer to additional data to be passed to the custom_y function
    void *custom_data_y;
    // Total number of particles already injected along x
    unsigned long custom_x_total_part;
    // Total charge injected (density integral) along x	
    double custom_x_total_q;

} t_density;


typedef struct {

    char name[MAX_SPNAME_LEN];

    // Particle data buffer
    t_part *part;
    int np;
    int np_max;

    // mass over charge ratio
    t_part_data m_q;

    // total kinetic energy
    double energy;

    // charge of individual particle
    t_part_data q;

    // Number of particles per cell
    int ppc[2];

    // Density profile to inject
    t_density density;

    // Initial momentum of particles
    t_part_data ufl[3];
    t_part_data uth[3];

    // Simulation box info
    int nx[2];
    t_part_data dx[2];
    t_part_data box[2];

    // Time step
    float dt;

    // Iteration number
    int iter;

    // Moving window
    int moving_window;
    int n_move;

} t_species;

void spec_new( t_species* spec, char name[], const t_part_data m_q, const int ppc[],
              const t_part_data ufl[], const t_part_data uth[],
              const int nx[], t_part_data box[], const float dt, t_density* density );

void spec_move_window( t_species *spec );

void spec_delete( t_species* spec );

void spec_advance( t_species* spec, t_emf* emf, t_current* current );

double spec_time( void );
double spec_perf( void );

/*********************************************************************************************

 Diagnostics

 *********************************************************************************************/

#define CHARGE 		0x1000
#define PHA    		0x2000
#define PARTICLES   0x3000
#define X1     		0x0001
#define X2     		0x0002
#define U1     		0x0004
#define U2     		0x0005
#define U3     		0x0006

#define PHASESPACE(a,b) ((a) + (b)*16 + PHA)

void spec_deposit_pha( const t_species *spec, const int rep_type,
              const int pha_nx[], const float pha_range[][2], float* buf );

void spec_report( const t_species *spec, const int rep_type,
                  const int pha_nx[], const float pha_range[][2] );

void spec_deposit_charge( const t_species* spec, float* charge );


#endif
