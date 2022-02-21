/*
 *  charge.h
 *  zpic
 *
 *  Created by Ricardo Fonseca on 12/8/10.
 *  Copyright 2010 Centro de Física dos Plasmas. All rights reserved.
 *
 */

#ifndef __CHARGE__
#define __CHARGE__

#include "zpic.h"
#include "grid.h"
#include "fft.h"
#include "filter.h"

/**
 * @brief Charge density object
 * 
 */
typedef struct Charge {
	
	/// Global charge density
	t_scalar_grid rho;

	/// Neutralizing background
	t_scalar_grid neutral;
	
	/// Fourier transform of rho
	t_cscalar_grid frho;

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

} t_charge;

/**
 * @brief Initializes Electric charge density object
 * 
 * @param charge 		Electric charge density
 * @param nx 			Number of cells
 * @param box 			Physical box size
 * @param dt 			Simulation time step
 * @param fft_forward 	FFT configuration for transforming rho to frho
 * 						(shared with other objects)
 * @param filter 		Spectral filtering parameters
 */
void charge_new( t_charge *charge, int nx, float box, float dt, t_fftr_cfg *fft_forward,
                 t_filter *filter );

/**
 * @brief Frees dynamic memory from electric charge density
 * 
 * @param charge 	Electric charge density
 */
void charge_delete( t_charge *charge );

/**
 * @brief Sets all electric charge density values to zero
 * 
 * @param charge 	Electric charge density
 */
void charge_zero( t_charge *charge );

/**
 * @brief Advances electric charge density 1 time step
 * 
 * @param charge 	Electric charge density
 */
void charge_update( t_charge *charge );

/**
 * @brief Saves electric charge density diagnostic information to disk
 * 
 * @param charge 	Electric charge density
 */
void charge_report( const t_charge *charge );

/**
 * @brief Initializes neutralizing background structures
 * 
 * @param charge 	Electric charge density
 */
void charge_init_neutral_bkg( t_charge *charge );

/**
 * @brief Updates neutralizing background values
 * 
 * @param charge 	Electric charge density
 */
void charge_update_neutral_bkg( t_charge *charge );


#endif
