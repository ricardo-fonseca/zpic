/*
 *  field.h
 *  zpic
 *
 *  Created by Ricardo Fonseca on 10/8/10.
 *  Copyright 2010 Centro de Física dos Plasmas. All rights reserved.
 *
 */

#ifndef __field__
#define __field__


#include "zpic.h"
#include "grid.h"
#include "charge.h"

/**
 * @brief Electric field
 * 
 */
typedef struct Field {
	
	/// E field
	t_scalar_grid E;

	/// Fourier transform of E
	t_cscalar_grid fE;
		
	// Simulation box info
	float box;	///< Physical size of simulation box
	float dx;	///< Grid cell size

	/// Time step
	float dt;

	/// Current Iteration number
	int iter;

	/// FFT configuration
	t_fftr_cfg *fft_backward;
	
} t_field;
	

/**
 * @brief Initializes the electric field object
 * 
 * @param field 			Electric field object
 * @param nx 				Number of cells
 * @param box 				Physical box size
 * @param dt 				Simulation time step
 * @param fft_backward 		FFT configuration (shared with other objects)
 */
void field_new( t_field *field, int nx, float box, const float dt, t_fftr_cfg *fft_backward );

/**
 * @brief Frees dynamic memory from electric field
 * 
 * @param field 	Electric field
 */
void field_delete( t_field *field );

/**
 * @brief Saves electric field diagnostic information to disk
 * 
 * Field will be save in directory "field", guard cell values are
 * discarded.
 * 
 * @param field 	Electric field
 */
void field_report( const t_field *field  );

/**
 * @brief Advance electric field 1 timestep
 * 
 * Field is updated from the charge density. The routine will also update
 * guard cell values.
 * 
 * @param field 
 * @param charge 
 */
void field_advance( t_field *field, const t_charge *charge );

/**
 * @brief Time spent advancing the electric fields
 * 
 * @return      Time spent in seconds
 */
double field_time( void );

#endif
