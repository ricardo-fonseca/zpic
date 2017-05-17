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

typedef struct {
	
	t_vfld_grid J;
	
	// Fourier transform of Jt
	t_cvfld_grid fJt;

	// Box size
	t_fld box;
	
	// Cell size
	t_fld dx;

	// Time step
	float dt;

	// Iteration number
	int iter;

	// FFT configuration
	t_fftr_cfg *fft_forward;
	
} t_current;

void current_new( t_current *current, int nx, t_fld box, float dt, t_fftr_cfg *fft );
void current_delete( t_current *current );
void current_zero( t_current *current );
void current_update( t_current *current );
void current_report( const t_current *current, const char jc );

#endif
