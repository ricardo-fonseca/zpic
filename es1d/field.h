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


typedef struct {
	
	// E field
	t_scalar_grid E;

	// Fourier transform of E
	t_cscalar_grid fE;
		
	// Simulation box info
	t_fld box;
	t_fld dx;

	// Time step
	float dt;

	// Iteration number
	int iter;

	// FFT configuration
	t_fftr_cfg *fft_backward;
	
} t_field;
	

void field_new( t_field *field, int nx, t_fld box, const float dt, t_fftr_cfg *fft_backward );

void field_delete( t_field *field );

void field_report( const t_field *field  );

void field_advance( t_field *field, const t_charge *charge );

void field_update_gc( t_field *field );

double field_time( void );

#endif
