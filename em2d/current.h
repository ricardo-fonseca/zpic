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

enum smooth_type { none, binomial, compensated };

typedef struct {
	enum smooth_type xtype, ytype;
	int xlevel, ylevel;
} t_smooth;

typedef struct {
	
	t_vfld *J;
	
	t_vfld *J_buf;
	
	// Grid parameters
	int nx[2];
	int nrow;
	int gc[2][2];
	
	// Box size
	t_fld box[2];
	
	// Cell size
	t_fld dx[2];

	// Current smoothing
	t_smooth smooth;

	// Time step
	float dt;

	// Iteration number
	int iter;

	// Moving window
	int moving_window;
	
} t_current;

void current_new( t_current *current, int nx[], t_fld box[], float dt );
void current_delete( t_current *current );
void current_zero( t_current *current );
void current_update( t_current *current );
void current_report( const t_current *current, const char jc );
void current_smooth( t_current* const current );

#endif
