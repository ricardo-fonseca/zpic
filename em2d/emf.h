/*
 *  emf.h
 *  zpic
 *
 *  Created by Ricardo Fonseca on 10/8/10.
 *  Copyright 2010 Centro de Física dos Plasmas. All rights reserved.
 *
 */

#ifndef __EMF__
#define __EMF__


#include "zpic.h"

#include "current.h"

enum emf_diag { EFLD, BFLD };

typedef struct {
	
	t_vfld *E;
	t_vfld *B;
	
	t_vfld *E_buf;
	t_vfld *B_buf;
	
	// Simulation box info
	int nx[2];
	int nrow;
	int gc[2][2];
	t_fld box[2];
	t_fld dx[2];

	// Time step
	float dt;

	// Iteration number
	int iter;

	// Moving window
	int moving_window;
	int n_move;
	
} t_emf;

enum emf_laser_type{ PLANE, GAUSSIAN };

typedef struct {

	enum emf_laser_type type;		// Laser pulse type
	
	float start;	// Front edge of the laser pulse, in simulation units
	float fwhm;		// FWHM of the laser pulse duration, in simulation units
	
	float a0;		// Normalized peak vector potential of the pulse
	float omega0;	// Laser frequency, normalized to the plasma frequency
	
	float polarization; 
	
	float W0;		// Gaussian beam waist, in simulation units
	float focus;	// Focal plane position, in simulation units
	float axis;     // Position of optical axis, in simulation units
	
} t_emf_laser;
	

void emf_new( t_emf *emf, int nx[], t_fld box[], const float dt );
void emf_delete( t_emf *emf );
void emf_report( const t_emf *emf, const char field, const char fc );

void emf_add_laser( t_emf* const emf, const t_emf_laser* const laser );

void emf_advance( t_emf *emf, const t_current *current );

void emf_move_window( t_emf *emf );

void emf_update_gc( t_emf *emf );

double emf_time();

void emf_set_moving_window( t_emf* emf );

#endif
