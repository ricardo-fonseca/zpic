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
#include "grid2d.h"
#include "fft.h"

/**
 * @brief Charge density
 * 
 */
typedef struct Charge {
	
	/// Global charge density
	t_scalar_grid2d rho;

	/// Neutralizing background
	t_scalar_grid2d neutral;
	
	/// Fourier transform of rho
	t_cscalar_grid2d frho;

	/// Box size [x,y]
	float box[2];
	
	/// Cell size
	float dx[2];

	/// Time step
	float dt;

	/// Iteration number
	int iter;

	/// FFT configuration
	t_fftr2d_cfg fft_forward;
	
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
 * @param filter 		Spectral filtering parameterso
 */
void charge_new( t_charge *charge, const int nx[], float box[], float dt);

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
 * The routine will:
 * 1. Update the guard cells
 * 2. Add neutralizing background (if configured)
 * 3. Get the fourier transform of the charge
 * 4. Apply spectral filtering
 * 
 * @param charge 	Electric charge density
 */
void charge_update( t_charge *charge );

/**
 * @brief Saves electric charge density diagnostic information to disk
 * 
 * Saves the charge density to disk in directory "CHARGE"
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
 * The routine expects the charge -> neutral grid to have the charge
 * density that needs to be neutralized. The routine will update guard
 * cell values and reverse the sign of the charge.
 * 
 * @param charge 
 */
void charge_update_neutral_bkg( t_charge *charge );


#endif
