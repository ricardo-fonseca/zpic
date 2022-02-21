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
#include "grid.h"
#include "fft.h"
#include "filter.h"

/**
 * @brief Current density object
 * 
 */
typedef struct Current {

	/// Current density grid
	t_float3_grid J;

	/// Fourier transform of J
	t_cfloat3_grid fJ;

	/// Box size
	float box;

	/// Cell size
	float dx;

	/// Time step
	float dt;

	/// Iteration number
	int iter;

	/// FFT configuration
	t_fftr_cfg *fft_forward;

	/// Spectral filtering
	t_filter *filter;

} t_current;

/**
 * @brief Initializes Electric current density object
 * 
 * @param current 		Electric current density
 * @param nx 			Number of cells
 * @param box 			Physical box size
 * @param dt 			Simulation time step
 * @param fft_forward 	FFT configuration for transforming J to fJ
 * 						(shared with other objects)
 * @param filter 		Spectral filtering parameters
 */
void current_new( t_current *current, int nx, float box, float dt, t_fftr_cfg *fft_forward,
	              t_filter *filter );

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
 * @param current Electric current density object
 */
void current_update( t_current *current );

/**
 * @brief Saves electric current density diagnostic information to disk
 * 
 * @param current Electric current density object
 * @param jc Current component to save, must be one of {0,1,2}
 */
void current_report( const t_current *current, const int jc );

#endif
