/*
 *  emf.h
 *  zpic
 *
 *  Created by Ricardo Fonseca on 10/8/10.
 *  Copyright 2010 Centro de Física dos Plasmas. All rights reserved.
 *
 */

#ifndef __EMF__
#define __EMF__

#include "zpic.h"
#include "current.h"

/**
 * @brief External/initial EM field types
 * 
 */
enum emf_fld_type { 
    EMF_FLD_TYPE_NONE,      ///< None
    EMF_FLD_TYPE_UNIFORM,   ///< Uniform
    EMF_FLD_TYPE_CUSTOM     ///< Defined from an external function
};

/**
 * @brief EM external field parameters
 * 
 */
typedef struct EMF_ExternalField {

    enum emf_fld_type E_type;   ///< Type of external E field
    enum emf_fld_type B_type;   ///< Type of external B field
    
    float3 E_0; ///< Value of uniform external E field
    float3 B_0; ///< Value of uniform external B field

    float3 (*E_custom)(int, float, void*);  ///< Custom external E-field function
    float3 (*B_custom)(int, float, void*);  ///< Custom external B-field function

    void *E_custom_data; ///< Additional data to be passed to the E_custom function
    void *B_custom_data; ///< Additional data to be passed to the B_custom function

    float3 *E_part_buf; ///< E field seen by particles
    float3 *B_part_buf; ///< B field seen by particles
} t_emf_ext_fld;

/**
 * @brief EM initial field parameters
 * 
 */
typedef struct EMF_InitialField {

    enum emf_fld_type E_type;   ///< Type of initial E field
    enum emf_fld_type B_type;   ///< Type of initial B field
    
    float3 E_0; ///< Value for uniform initial E field
    float3 B_0; ///< Value for uniform initial B field

    float3 (*E_custom)(int, float, void*); ///< Custom initial E-field function
    float3 (*B_custom)(int, float, void*); ///< Custom initial B-field function

    void *E_custom_data; ///< Additional data to be passed to the E_custom function
    void *B_custom_data; ///< Additional data to be passed to the B_custom function

} t_emf_init_fld;

/**
 * @brief Types of EM field diagnostics
 * 
 */
enum emf_diag {
    EFLD,   ///< Electric field
    BFLD,   ///< Magnetic field
    EPART,  ///< Electric field as seen by the particles
    BPART   ///< Magnetic field as seen by the particles
};

/**
 * @brief Types of EM field boundaries
 * 
 */
enum emf_boundary {
    EMF_BC_NONE,        ///< None
    EMF_BC_PERIODIC,    ///< Periodic boundaries
    EMF_BC_OPEN         ///< Open (Mur type) boundaries
};


/**
 * @brief Electro-Magnetic fields
 * 
 */
typedef struct EMF {

    float3 *E;  ///< Pointer to grid cell 0 of E field
    float3 *B;  ///< Pointer to grid cell 0 of B field

    float3 *E_buf;  ///< E field buffer (includes guard cells)
    float3 *B_buf;  ///< B field buffer (includes guard cells)

    // Fields seen by particles
    // When using external fields these will be a combination of the simulation
    // fields and the externally imposed ones. When external fields are off
    // these just point to E and B.
    
    float3 *E_part; ///< Pointer to grid cell 0 of particles E field
    float3 *B_part; ///< Pointer to grid cell 0 of particles B field

    // Simulation box info
    int nx;  ///< Number of grid points (excluding guard cells) 
    int gc[2];  ///< Number of guard cells (lower/upper)
    float box;  ///< Physical size of simulation box
    float dx;   ///< Grid cell size

    float dt;   ///< Time step

    int iter;   ///< Current iteration number

    // Moving window
    int moving_window;  ///< Active moving window
    int n_move;         ///< Number of cells moved by the moving window algorithm

    // Boundary conditions
    enum emf_boundary bc_type;  ///< Type of boundary conditions
    float3 mur_fld[2]; ///< Field inside Mur boundary (lower/upper)
    float3 mur_tmp[2]; ///< Previous Eperp for Mur boundary advance (lower/upper)

    /// External fields configuration
    t_emf_ext_fld ext_fld;

} t_emf;


/**
 * @brief Laser Pulse parameters
 * 
 */
typedef struct EMF_Laser {
    
    float start; ///< Front edge of the laser pulse, in simulation units
    
    float fwhm; ///< FWHM of the laser pulse duration, in simulation units

    float rise; ///< Rise time of the laser pulse, in simulation units
    float flat; ///< Flat time of the laser pulse, in simulation units
    float fall; ///< Fall time of the laser pulse, in simulation units

    float a0;		///< Normalized peak vector potential of the pulse
    float omega0;	///< Laser frequency, normalized to the simulation reference frequency

    float polarization; ///< Polarization angle in radians

} t_emf_laser;

/**
 * @brief Calculate total EM field energy
 *
 * @param[in] emf EM field
 * @param[out] energy Energy values vector
 */
void emf_get_energy( const t_emf *emf, double energy[] );

/**
 * @brief Initalized EM fields object
 * 
 * @param emf 	EM fields
 * @param nx 	Number of grid cells
 * @param box 	Physical box size
 * @param dt 	Simulation time step
 */
void emf_new( t_emf *emf, int nx, float box, const float dt );

/**
 * @brief Frees dynamic memory from EM fields.
 * 
 * If external fields are in use, the dynamic memory associated with these
 * will also be freed.
 * 
 * @param emf 	EM fields
 */
void emf_delete( t_emf *emf );

/**
 * @brief Saves EM fields diagnostic information to disk
 * 
 * Saves the selected type / density component to disk in directory
 * "EMF". Guard cell values are discarded.
 *
 * @param emf 		EM Fields
 * @param field 	Which field to save (E, B, Epart, Bpart)
 * @param fc 		Field component to save, must be one of {0,1,2}
 */
void emf_report( const t_emf *emf, const char field, const int fc );

/**
 * @brief Add laser pulse to simulation.
 * 
 * Laser pulses are superimposed on top of existing E and B fields. 
 * Multiple lasers can be added.
 * 
 * @param emf 		EM fields
 * @param laser 	Laser pulse parameters
 */
void emf_add_laser( t_emf* const emf, t_emf_laser* laser );

/**
 * @brief Initialize EMF field values
 * 
 * @param emf       EM field object
 * @param inifloat  Initial field parameters
 */
void emf_init_fld( t_emf* const emf, t_emf_init_fld* inifloat );

/**
 * @brief Sets the external fields to be used for the simulation
 * 
 * @param emf 		EM field
 * @param exfloat 	External fields
 */
void emf_set_ext_fld( t_emf* const emf, t_emf_ext_fld* exfloat );

/**
 * @brief Advance EM fields 1 timestep
 * 
 * Fields are advanced in time using a FDTD algorith. The routine will also:
 * 1. Update guard cell values / apply boundary conditions
 * 2. Update "particle" fields if using external fields
 * 3. Move simulation window 
 * 
 * @param emf 		EM fields
 * @param current 	Electric current density
 */
void emf_advance( t_emf *emf, const t_current *current );

/**
 * @brief Time spent advancing the EM fields
 * 
 * @return double 	Time spent in seconds
 */
double emf_time( void );

#endif
