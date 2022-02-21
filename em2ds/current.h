/*
 *  current.h
 *  zpic
 *
 *  Created by Ricardo Fonseca on 12/8/10.
 *  Copyright 2010 Centro de Física dos Plasmas. All rights reserved.
 *
 */

#ifndef __CURRENT__
#define __CURRENT__

#include "zpic.h"
#include "grid2d.h"
#include "fft.h"

/**
 * @brief Current density object
 * 
 */
typedef struct Current {

	/// Current density grid
	t_float3_grid2d J;

	/// Fourier transform of J
	t_cfloat3_grid2d fJ;

	/// Simulation box size [x,y]
	float box[2];

	/// Cell size [x,y]
	float dx[2];

	/// Time step
	float dt;

	/// Current iteration number
	int iter;

	// FFT configuration
	t_fftr2d_cfg fft_forward;	///< Forward (real to complex) 2D FFT configuration
	t_fftr2d_cfg fft_backward;	///< Backward (complex to real) 2D FFT configuration

} t_current;

/**
 * @brief Initializes Electric current density object
 * 
 * @param current 		Electric current density
 * @param nx 			Number of cells
 * @param box 			Physical box size
 * @param dt 			Simulation time step
 */
void current_new( t_current *current, const int nx[], float box[], float dt);

/**
 * @brief Frees dynamic memory from electric current density
 * 
 * @param current Electric current density object
 */
void current_delete( t_current *current );

/**
 * @brief Sets all electric current density values to zero
 * 
 * @param current Electric current density object
 */
void current_zero( t_current *current );

/**
 * @brief Advances electric current density 1 time step
 * 
 * The routine will:
 * 1. Update the guard cells
 * 2. Calculate the transverse component of J (not needed in 1D)
 * 3. Get the Fourier transform of the current
 * 4. Apply spectral filtering
 * 
 * @param current Electric current density object
 */
void current_update( t_current *current );

/**
 * @brief Saves electric current density diagnostic information to disk
 * 
 * Saves the selected current density component to disk in directory
 * "CURRENT". Guard cell values are discarded.
 * 
 * @param current Electric current object
 * @param jc Current component to save, must be one of {0,1,2}
 */
void current_report( const t_current *current, const char jc );

#endif
