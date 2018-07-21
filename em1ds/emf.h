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
#include "grid.h"
#include "current.h"
#include "charge.h"

enum emf_ext_fld { EMF_EXT_FLD_NONE, EMF_EXT_FLD_UNIFORM };

/*
Some headers (e.g. termios.h) define a B0 macro that conflicts with this code
 */
#ifdef B0
#undef B0
#endif

typedef struct {
	t_vfld E0;
	t_vfld B0;
	enum emf_ext_fld type;

	t_vfld_grid E_part_buf;
	t_vfld_grid B_part_buf;
} t_emf_ext_fld;

enum emf_diag { EFLD, BFLD };

enum emf_solver { EMF_SOLVER_PSTD, EMF_SOLVER_PSATD };

typedef struct {

	// E and B fields
	t_vfld_grid E;
	t_vfld_grid B;

	// Fourier transform of El, Et and B
	t_cvfld_grid fEl;
	t_cvfld_grid fEt;
	t_cvfld_grid fB;

	// Simulation box info
	t_fld box;
	t_fld dx;

	// Time step
	float dt;

	// Iteration number
	int iter;

	// FFT configuration
	t_fftr_cfg *fft_forward, *fft_backward;

	// Fields seen by particles
	// When using external fields these will be a combination of the simulation
	// fields and the externally imposed ones. When external fields are off
	// these just point to E and B.
	t_vfld_grid *E_part;
	t_vfld_grid *B_part;

	// External fields
	t_emf_ext_fld ext_fld;

	// solver type
	enum emf_solver solver_type;

} t_emf;


// In 1D we only have plane waves
typedef struct {

	float start;	// Front edge of the laser pulse, in simulation units
	float fwhm;		// FWHM of the laser pulse duration, in simulation units
	float rise, flat, fall; // Rise, flat and fall time of the laser pulse, in simulation units

	float a0;		// Normalized peak vector potential of the pulse
	float omega0;	// Laser frequency, normalized to the plasma frequency

	float polarization;

} t_emf_laser;


void emf_new( t_emf *emf, int nx, t_fld box, const float dt, t_fftr_cfg *fft_forward,
	t_fftr_cfg *fft_backward );

void emf_delete( t_emf *emf );

void emf_report( const t_emf *emf, const char field, const char fc );

void emf_add_laser( t_emf* const emf, const t_emf_laser* const laser );

void emf_advance( t_emf *emf, const t_charge *charge, const t_current *current );

void emf_update_part_fld( t_emf *emf );

void emf_move_window( t_emf *emf );

void emf_update_gc( t_emf *emf );

double emf_time( void );

void emf_get_energy( const t_emf *emf, double energy[] );

void emf_set_ext_fld( t_emf* const emf, t_emf_ext_fld* ext_fld );


#endif
