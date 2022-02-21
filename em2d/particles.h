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
    char name[MAX_SPNAME_LEN];

    // Particles
    t_part *part;   ///< Particle buffer
    int np;         ///< Number of particles in buffer
    int np_max;     ///< Maximum number of particles in buffer

    /// mass over charge ratio
    float m_q;

    /// total kinetic energy
    double energy;

    /// charge of individual particle
    float q;

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

	// Moving window
	int moving_window;  ///< Active moving window
	int n_move;         ///< Number of cells moved by the moving window algorithm

	/// Sorting frequency
	int n_sort;

} t_species;

/**
 * @brief Initialize particle Species object
 * 
 * This routine will also inject the initial particle distribution,
 * setting thermal/fluid velocities.
 *  
 * @param spec      Particle species
 * @param name      Name for the species (used for diagnostic output)
 * @param m_q       Mass over charge ratio for species, in simulation units
 * @param ppc       Reference number of particles per cell [x,y]
 * @param ufl       Initial fluid momentum of particles, may be set to NULL 
 * @param uth       Initial thermal momentum of particles, may be set to NULL
 * @param nx        Number of grid points [x,y]
 * @param box       Simulation box size [x,y] in simulation units
 * @param dt        Simulation time step, in simulation units
 * @param density   Density profile for particle injection, may be set to NULL
 */
void spec_new( t_species* spec, char name[], const float m_q, const int ppc[],
              const float ufl[], const float uth[],
              const int nx[], float box[], const float dt, t_density* density );

/**
 * @brief Frees dynamic memory from particle species
 * 
 * @param spec Particle species
 */
void spec_delete( t_species* spec );

/**
 * @brief Grows particle buffer to specified size.
 * 
 * If the new size is smaller than the previous size the buffer size is not changed
 * and the function returns silently.
 * 
 * @param   spec    Particle species
 * @param   size    New buffer size (will be rounded up to next multiple of 1024)
 */
void spec_grow_buffer( t_species* spec, const int size );

/**
 * @brief Advance Particle species 1 timestep
 * 
 * Particles are advanced in time using a leap-frog method; the velocity
 * advance is done using a relativistic Boris pusher. Particle motion is
 * used to deposit electric current on the grid using an exact charge
 * conservation method.
 * 
 * The routine will also:
 * 1. Calculate total time-centered kinetic energy for the Species
 * 2. Apply boundary conditions
 * 3. Move simulation window
 * 4. Sort particle buffer
 * 
 * @param spec      Particle species
 * @param emf       EM fields
 * @param current   Current density
 */
void spec_advance( t_species* spec, t_emf* emf, t_current* current );

/**
 * @brief Move simulation window
 * 
 * When using a moving simulation window checks if a window move is due
 * at the current iteration and if so shifts left the particle cell indices
 * 
 * @param spec      Particle species
 */
void spec_move_window( t_species *spec );

/**
 * @brief Returns the total number of particle pushes
 * 
 * @return  Number of particle pushes
 */
uint64_t spec_npush( void );

/**
 * @brief Returns the total time spent pushing particles (includes boundaries and moving window)
 * 
 * @return  Total time in seconds
 */
double spec_time( void );

/**
 * @brief Returns the performance achieved by the code (push time)
 * 
 * @return  Performance in seconds per particle
 */
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

/**
 * @brief Defines a phasespace density report with the specified quantities
 * 
 * @param a		Phasespace x axis
 * @param b		Phasespace y axis
 * @return int value describing phasespace report
 */
#define PHASESPACE(a,b) ((a) + (b)*16 + PHA)

/**
 * @brief Deposit 2D phasespace density.
 * 
 * @param spec      Particle species
 * @param rep_type  Type of phasespace, use the PHASESPACE macro to define
 * @param pha_nx    Number of grid points in the phasespace density grid
 * @param pha_range Physical range of each of the phasespace axis
 * @param buf       Phasespace density grid
 */
void spec_deposit_pha( const t_species *spec, const int rep_type,
              const int pha_nx[], const float pha_range[][2], float* buf );

/**
 * @brief Saves particle species diagnostic information to disk
 * 
 * @param spec      Particle species
 * @param rep_type  Type of diagnostic information {CHARGE, PHASESPACE(a,b), PARTICLES}
 * @param pha_nx    Number of grid points in the phasespace density grid, set to NULL for
 *                  diagnostics other then phasespace density
 * @param pha_range Physical range of each of the phasespace axis, set to NULL for
 *                  diagnostics other then phasespace density
 */
void spec_report( const t_species *spec, const int rep_type,
                  const int pha_nx[], const float pha_range[][2] );

/**
 * @brief Deposits particle species charge density
 * 
 * Deposition is done using linear interpolation. Used for diagnostics
 * purpose only.
 * 
 * @param spec      Particle species
 * @param charge    Electric charge density
 */
void spec_deposit_charge( const t_species* spec, float* charge );


#endif
