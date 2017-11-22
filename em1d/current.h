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
enum current_boundary{ CURRENT_BC_NONE, CURRENT_BC_PERIODIC };

typedef struct {
	enum smooth_type xtype;
	int xlevel;
} t_smooth;

typedef struct {
	
	t_vfld *J;
	
	t_vfld *J_buf;
	
	// Grid parameters
	int nx;
	int gc[2];
	
	// Box size
	t_fld box;
	
	// Cell size
	t_fld dx;

	// Current smoothing
	t_smooth smooth;

	// Time step
	float dt;

	// Iteration number
	int iter;

	// Boundary conditions
	int bc_type;
	
} t_current;

void current_new( t_current *current, int nx, t_fld box, float dt );
void current_delete( t_current *current );
void current_zero( t_current *current );
void current_update( t_current *current );
void current_report( const t_current *current, const char jc );
void current_smooth( t_current* const current );

#endif
