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

#include <complex.h>
#include "zpic.h"

typedef struct {

	float *s;

	float *buffer;

	int nx[2];
	int nrow;
	int gc[2][2];

} t_scalar_grid2d;

typedef struct {

	float complex *s;

	float complex *buffer;

	int nx[2];
	int nrow;
	int gc[2][2];

} t_cscalar_grid2d;


typedef struct {
	float x, y, z;
} t_vfld;

typedef struct {

	float *x,*y,*z;

	float *buffer;

	int nx[2];
	int nrow;
	int gc[2][2];

} t_vfld_grid2d;


typedef struct {

	float complex *x,*y,*z;

	float complex *buffer;

	int nx[2];
	int nrow;
	int gc[2][2];

} t_cvfld_grid2d;


int  scalar_grid2d_init( t_scalar_grid2d *grid, const unsigned int nx[], const unsigned int * gc );
int  scalar_grid2d_cleanup( t_scalar_grid2d *grid );
void scalar_grid2d_zero( t_scalar_grid2d *grid );
void scalar_grid2d_copy( t_scalar_grid2d *dst, const t_scalar_grid2d *src );

int  cscalar_grid2d_init( t_cscalar_grid2d *grid, const unsigned int nx[], const unsigned int * gc );
int  cscalar_grid2d_cleanup( t_cscalar_grid2d *grid );
void cscalar_grid2d_zero( t_cscalar_grid2d *grid );


int  vfld_grid2d_init( t_vfld_grid2d *grid, const unsigned int nx[], const unsigned int * gc );
int  vfld_grid2d_cleanup( t_vfld_grid2d *grid );
void vfld_grid2d_zero( t_vfld_grid2d *grid );

int  cvfld_grid2d_init( t_cvfld_grid2d *grid, const unsigned int nx[], const unsigned int * gc );
int  cvfld_grid2d_cleanup( t_cvfld_grid2d *grid );
void cvfld_grid2d_zero( t_cvfld_grid2d *grid );

#endif
