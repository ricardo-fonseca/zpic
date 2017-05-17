/*
 *  fld_grid.h
 *  zpic
 *
 *  Created by Ricardo Fonseca on 12/8/10.
 *  Copyright 2010 Centro de Física dos Plasmas. All rights reserved.
 *
 */

#ifndef __GRID__
#define __GRID__

#include <complex.h>
#include "zpic.h"

typedef struct {
	
	float *s;
	
	float *buffer;
	
	int nx;
	int gc[2];
		
} t_scalar_grid;

typedef struct {
	
	float complex *s;
	
	float complex *buffer;
	
	int nx;
	int gc[2];
		
} t_cscalar_grid;


typedef struct {
	float x, y, z;
} t_vfld;

typedef struct {
	
	float *x,*y,*z;
	
	float *buffer;
	
	int nx;
	int gc[2];
		
} t_vfld_grid;


typedef struct {
	
	float complex *x,*y,*z;
	
	float complex *buffer;
	
	int nx;
	int gc[2];
		
} t_cvfld_grid;


int  scalar_grid_init( t_scalar_grid *grid, const unsigned int nx, const unsigned int * gc );
int  scalar_grid_cleanup( t_scalar_grid *grid );
void scalar_grid_zero( t_scalar_grid *grid );
void scalar_grid_copy( t_scalar_grid *dst, t_scalar_grid *src );

int  cscalar_grid_init( t_cscalar_grid *grid, const unsigned int nx, const unsigned int * gc );
int  cscalar_grid_cleanup( t_cscalar_grid *grid );
void cscalar_grid_zero( t_cscalar_grid *grid );


int  vfld_grid_init( t_vfld_grid *grid, const unsigned int nx, const unsigned int * gc );
int  vfld_grid_cleanup( t_vfld_grid *grid );
void vfld_grid_zero( t_vfld_grid *grid );

int  cvfld_grid_init( t_cvfld_grid *grid, const unsigned int nx, const unsigned int * gc );
int  cvfld_grid_cleanup( t_cvfld_grid *grid );
void cvfld_grid_zero( t_cvfld_grid *grid );

#endif
