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
#include "charge.h"

#include <stdint.h>

/**
 * @brief Maximum species name length
 * 
 */
#define MAX_SPNAME_LEN 32

/**
 * @brief Structure holding single particle data
 * 
 */
typedef struct Particle {
	int ix;		///< Particle x cell index
	int iy;		///< Particle y cell index
	float x;	///< x position inside cell
	float y;	///< y position inside cell
	float ux;	///< Generalized velocity along x
	float uy;	///< Generalized velocity along y
	float uz;	///< Generalized velocity along z
} t_part;

/**
 * @brief Types of density profile
 * 
 */
enum density_type {
	UNIFORM,	///< Uniform density
	EMPTY,		///< No particles
	STEP,		///< Step-like profile
	SLAB,		///< Slab-like profile
	CUSTOM		///< Defined from an external function
};

/**
 * @brief Density profile parameters
 * 
 */
typedef struct Density {

	float n;						///< reference density (defaults to 1.0, multiplies density profile)

	enum density_type type;			///< Density profile type
	float start;					///< Start position for step, slab and ramp profiles, in simulation units
	float end;						///< End position for slab and ramp profiles, in simulation units

    // Custom density profile parameters
    
    /// Pointer to custom density function along x
    float (*custom_x)(float, void*);
    /// Pointer to additional data to be passed to the custom_x function
    void *custom_data_x;
    /// Pointer to custom density function along y
    float (*custom_y)(float, void*);
    /// Pointer to additional data to be passed to the custom_y function
    void *custom_data_y;
    /// Total number of particles already injected along x
    unsigned long custom_x_total_part;
    /// Total charge injected (density integral) along x
    double custom_x_total_q;

} t_density;

/**
 * @brief Set of particles
 * 
 */
typedef struct Species {
    
    /// Species name
    char name[MAX_SPNAME_LEN+1];
    
    // Particles
    t_part *part;   ///< Particle buffer
    int np;         ///< Number of particles in buffer
    int np_max;     ///< Maximum number of particles in buffer

    /// mass over charge ratio
    float m_q;
    
    /// charge of individual particle
    float q;

    /// total kinetic energy
    double energy;

    /// Number of particles per cell [x,y]
    int ppc[2];

    /// Density profile to inject
    t_density density;

    // Initial momentum of particles
    float ufl[3];   ///< Initial fluid momentum of particles
    float uth[3];   ///< Initial thermal momentum of particles

    // Simulation box info
    int nx[2];      ///< Number of grid points [x,y] (excluding guard cells)
    float dx[2];    ///< Grid cell size [x,y]
    float box[2];   ///< Physical size of simulation box [x,y]

    /// Time step
    float dt;

    /// Current iteration number
    int iter;

	/// Sorting frequency
	int n_sort;

} t_species;

void spec_new( t_species* spec, char name[], const float m_q, const int ppc[], 
              const float ufl[], const float uth[],
              const int nx[], float box[], const float dt, t_density* density );

void spec_delete( t_species* spec );

void spec_grow_buffer( t_species* spec, const int size );

void spec_advance( t_species* spec, t_emf* emf, t_charge* charge, t_current* current );

void spec_deposit_charge( const t_species* spec, float* charge );

uint64_t spec_npush( void );
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

void spec_deposit_charge( const t_species* spec, float* charge );

void spec_report( const t_species *spec, const int rep_type,
                  const int pha_nx[], const float pha_range[][2] );


#endif
