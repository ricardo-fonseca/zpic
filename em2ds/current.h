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


typedef struct {

	t_vfld_grid2d J;

	// Fourier transform of Jt
	t_cvfld_grid2d fJt;

	// Box size
	t_fld box[2];

	// Cell size
	t_fld dx[2];

	// Time step
	float dt;

	// Iteration number
	int iter;

	// FFT configuration
	t_fftr2d_cfg fft_forward, fft_backward;

} t_current;

void current_new( t_current *current, const int nx[], t_fld box[], float dt);
void current_delete( t_current *current );
void current_zero( t_current *current );
void current_update( t_current *current );
void current_report( const t_current *current, const char jc );

#endif
