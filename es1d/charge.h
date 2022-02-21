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

/**
 * @brief Charge density
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
	
} t_charge;


void charge_new( t_charge *charge, int nx, float box, float dt, t_fftr_cfg *fft_forward );
void charge_delete( t_charge *charge );
void charge_zero( t_charge *charge );
void charge_update( t_charge *charge );
void charge_report( const t_charge *charge );

void charge_init_neutral_bkg( t_charge *charge );
void charge_update_neutral_bkg( t_charge *charge );


#endif
