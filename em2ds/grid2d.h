/*
 *  fld_grid.h
 *  zpic
 *
 *  Created by Ricardo Fonseca on 12/8/10.
 *  Copyright 2010 Centro de Física dos Plasmas. All rights reserved.
 *
 */

#ifndef __GRID_2D__
#define __GRID_2D__

#include <math.h>
#include <complex.h>
#include "zpic.h"

typedef struct {

	float *s;

	float *buffer;

	unsigned int nx[2];
	unsigned int nrow;
	unsigned int gc[2][2];

} t_scalar_grid2d;

typedef struct {

	float complex *s;

	float complex *buffer;

	unsigned int nx[2];
	unsigned int nrow;
	unsigned int gc[2][2];

} t_cscalar_grid2d;

typedef struct {

	float *x,*y,*z;

	float *buffer;

	unsigned int nx[2];
	unsigned int nrow;
	unsigned int gc[2][2];

} t_vfld_grid2d;


typedef struct {

	float complex *x,*y,*z;

	float complex *buffer;

	unsigned int nx[2];
	unsigned int nrow;
	unsigned int gc[2][2];

} t_cvfld_grid2d;


int  scalar_grid2d_init( t_scalar_grid2d *grid, const unsigned int nx[], const unsigned int gc[][2] );
int  scalar_grid2d_cleanup( t_scalar_grid2d *grid );
void scalar_grid2d_zero( t_scalar_grid2d *grid );
void scalar_grid2d_copy( t_scalar_grid2d * restrict dst, const t_scalar_grid2d * restrict src );
void scalar_grid2d_add( t_scalar_grid2d * restrict dst, const t_scalar_grid2d * restrict src );


int cscalar_grid2d_init( t_cscalar_grid2d *grid, const unsigned int nx[], const unsigned int gc[][2] );
int  cscalar_grid2d_cleanup( t_cscalar_grid2d *grid );
void cscalar_grid2d_zero( t_cscalar_grid2d *grid );
void cscalar_grid2d_copy( t_cscalar_grid2d * restrict dst, const t_cscalar_grid2d * restrict src );
void cscalar_grid2d_add( t_cscalar_grid2d * restrict dst, const t_cscalar_grid2d * restrict src );


int  vfld_grid2d_init( t_vfld_grid2d *grid, const unsigned int nx[], const unsigned int gc[][2] );
int  vfld_grid2d_cleanup( t_vfld_grid2d *grid );
void vfld_grid2d_zero( t_vfld_grid2d *grid );
void vfld_grid2d_copy( t_vfld_grid2d * restrict dst, const t_vfld_grid2d * restrict src );
void vfld_grid2d_add( t_vfld_grid2d * restrict dst, const t_vfld_grid2d * restrict src );


int  cvfld_grid2d_init( t_cvfld_grid2d *grid, const unsigned int nx[], const unsigned int gc[][2] );
int  cvfld_grid2d_cleanup( t_cvfld_grid2d *grid );
void cvfld_grid2d_zero( t_cvfld_grid2d *grid );
void cvfld_grid2d_copy( t_cvfld_grid2d * restrict dst, const t_cvfld_grid2d * restrict src );
void cvfld_grid2d_add( t_cvfld_grid2d * restrict dst, const t_cvfld_grid2d * restrict src );

#endif
