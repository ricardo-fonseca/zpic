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
#include "grid2d.h"
#include "current.h"
#include "charge.h"

/**
 * @brief Types of EM field diagnostics
 * 
 */
enum emf_diag {
    EFLD,   ///< Electric field
    BFLD    ///< Magnetic field
};

/**
 * @brief Types of EM field solvers
 * 
 */
enum emf_solver {
	EMF_SOLVER_PSTD,	///< Pseudo-spectral time domain
	EMF_SOLVER_PSATD	///< Pseudo-spectral analytic time domain
};

/**
 * @brief EM Fields
 * 
 */
typedef struct EMF {

	// E and B fields
	t_float3_grid2d E;		///< Electric field grid
	t_float3_grid2d B;		///< Magnetic field grid

	t_cfloat3_grid2d fEl;	///< Fourier transform of longitudinal Electric field
	t_cfloat3_grid2d fEt;	///< Fourier transform of transverse Electric field
	t_cfloat3_grid2d fB;	///< Fourier transform of Magnetic field

	// Simulation box info
	float box[2];	///< Physical size of simulation box [x,y]
	float dx[2];	///< Grid cell size [x,y]

	float dt;		///< Time step
	int iter;		///< Current iteration number

	// FFT configurations
	t_fftr2d_cfg fft_forward;	///< FFT configuration for forward transformation
	t_fftr2d_cfg fft_backward;	///< FFT configuration for backward transformation

	/// EM solver type
	enum emf_solver solver_type;

} t_emf;

/**
 * @brief Types of laser pulses
 * 
 */
enum emf_laser_type{
	PLANE,		///< Plane wave
	GAUSSIAN	///< Gaussian beam
};

/**
 * @brief Laser Pulse parameters
 * 
 */
typedef struct EMF_Laser {

	enum emf_laser_type type;  ///< Laser pulse type
	
	float start;  ///< Front edge of the laser pulse, in simulation units

	float fwhm;  ///< FWHM of the laser pulse duration, in simulation units

    float rise; ///< Rise time of the laser pulse, in simulation units
    float flat; ///< Flat time of the laser pulse, in simulation units
    float fall; ///< Fall time of the laser pulse, in simulation units
	
	float a0;		///< Normalized peak vector potential of the pulse
	float omega0;	///< Laser frequency, normalized to the plasma frequency
	
	float polarization; ///< Polarization angle in radians
	
	float W0;		///< Gaussian beam waist, in simulation units
	float focus;	///< Focal plane position, in simulation units
	float axis;     ///< Position of optical axis, in simulation units
	
} t_emf_laser;

void emf_get_energy( const t_emf *emf, double energy[] );

void emf_new( t_emf *emf, const int nx[], float box[], const float dt );
void emf_delete( t_emf *emf );
void emf_report( const t_emf *emf, const char field, const int fc );

void emf_add_laser( t_emf* const emf, const t_emf_laser* const laser );

void emf_advance( t_emf *emf, const t_charge *charge, const t_current *current );

void emf_update_gc( t_emf *emf );

double emf_time( void );

#endif
