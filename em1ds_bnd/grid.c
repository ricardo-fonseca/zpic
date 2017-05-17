/*
 *  grid.c
 *  zpic
 *
 *  Created by Ricardo Fonseca on 12/8/10.
 *  Copyright 2010 Centro de Física dos Plasmas. All rights reserved.
 *
 */

#include "grid.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/****************************************************************************************
	Scalar grids
 ****************************************************************************************/

int scalar_grid_init( t_scalar_grid *grid, const unsigned int nx, const unsigned int * gc )
{

	// store nx and gc values
	grid->nx = nx;
	if ( gc ) {
		grid->gc[0] = gc[0];
		grid->gc[1] = gc[1];
	} else {
		grid->gc[0] = 0;
		grid->gc[1] = 0;
	}


	size_t size = grid->gc[0] + grid->nx + grid->gc[1];
	
	grid -> buffer = (float *) malloc( size * sizeof( float ) );
	
	if ( !grid -> buffer ) {
		fprintf(stderr, "(*error*) Unable to allocate memory for fld_grid variable\n");
		return(-1);
	}

	
	// Make s point to cell [0]
	grid -> s = grid -> buffer + grid->gc[0];

	return 0;
}


int scalar_grid_cleanup( t_scalar_grid *grid )
{
	free( grid -> buffer );

	grid -> buffer = NULL;
	grid -> s = NULL;

	grid->nx = 0;
	grid->gc[0] = 0;
	grid->gc[1] = 0;

	return 0;

}

void scalar_grid_zero( t_scalar_grid *grid ) {
	size_t size =  (grid->gc[0] + grid->nx + grid->gc[1]) * sizeof(float);
	memset( grid -> buffer, 0, size );
}

void scalar_grid_copy( t_scalar_grid *dst, t_scalar_grid *src ){
	size_t n = ( dst -> gc[0] + dst -> nx + dst -> gc[1] ) * sizeof( float );
	memcpy( dst->buffer, src->buffer, n);
}


/****************************************************************************************
	Complex Scalar grids
 ****************************************************************************************/

int cscalar_grid_init( t_cscalar_grid *grid, const unsigned int nx, const unsigned int * gc )
{

	// store nx and gc values
	grid->nx = nx;
	if ( gc ) {
		grid->gc[0] = gc[0];
		grid->gc[1] = gc[1];
	} else {
		grid->gc[0] = 0;
		grid->gc[1] = 0;
	}


	size_t size = grid->gc[0] + grid->nx + grid->gc[1];
	
	grid -> buffer = (float complex *) malloc( size * sizeof( float complex ) );
	
	if ( !grid -> buffer ) {
		fprintf(stderr, "(*error*) Unable to allocate memory for fld_grid variable\n");
		return(-1);
	}

	
	// Make s point to cell [0]
	grid -> s = grid -> buffer + grid->gc[0];

	return 0;
}


int cscalar_grid_cleanup( t_cscalar_grid *grid )
{
	free( grid -> buffer );

	grid -> buffer = NULL;
	grid -> s = NULL;

	grid->nx = 0;
	grid->gc[0] = 0;
	grid->gc[1] = 0;

	return 0;

}

void cscalar_grid_zero( t_cscalar_grid *grid ) {
	size_t size =  (grid->gc[0] + grid->nx + grid->gc[1]) * sizeof(float complex);
	memset( grid -> buffer, 0, size );
}

/****************************************************************************************
	Vector field grids
 ****************************************************************************************/

int vfld_grid_init( t_vfld_grid *grid, const unsigned int nx, const unsigned int * gc )
{

	// store nx and gc values
	grid->nx = nx;
	if ( gc ) {
		grid->gc[0] = gc[0];
		grid->gc[1] = gc[1];
	} else {
		grid->gc[0] = 0;
		grid->gc[1] = 0;
	}


	size_t size = grid->gc[0] + grid->nx + grid->gc[1];
	
	grid -> buffer = malloc( 3 * size * sizeof( float ) );
	
	if ( !grid -> buffer ) {
		fprintf(stderr, "(*error*) Unable to allocate memory for fld_grid variable\n");
		return(-1);
	}

	
	// Make x, y and z point to cell [0]
	grid -> x = grid -> buffer + grid->gc[0];
	grid -> y = grid -> x + size;
	grid -> z = grid -> y + size;

	return 0;
}


int vfld_grid_cleanup( t_vfld_grid *grid )
{
	free( grid -> buffer );

	grid -> buffer = NULL;
	grid -> x = grid -> y = grid -> z = NULL;
	grid->nx = 0;
	grid->gc[0] = 0;
	grid->gc[1] = 0;

	return 0;

}

void vfld_grid_zero( t_vfld_grid *grid ) {
	size_t size = 3 * (grid->gc[0] + grid->nx + grid->gc[1]) * sizeof(float);
	memset( grid -> buffer, 0, size );
}

/****************************************************************************************
	Complex Vector field grids
 ****************************************************************************************/

int cvfld_grid_init( t_cvfld_grid *grid, const unsigned int nx, const unsigned int * gc )
{

	// store nx and gc values
	grid->nx = nx;
	if ( gc ) {
		grid->gc[0] = gc[0];
		grid->gc[1] = gc[1];
	} else {
		grid->gc[0] = 0;
		grid->gc[1] = 0;
	}


	size_t size = grid->gc[0] + grid->nx + grid->gc[1];
	
	grid -> buffer = malloc( 3 * size * sizeof( float complex ) );
	
	if ( !grid -> buffer ) {
		fprintf(stderr, "(*error*) Unable to allocate memory for fld_grid variable\n");
		return(-1);
	}

	
	// Make x, y and z point to cell [0]
	grid -> x = grid -> buffer + grid->gc[0];
	grid -> y = grid -> x + size;
	grid -> z = grid -> y + size;

	return 0;
}


int cvfld_grid_cleanup( t_cvfld_grid *grid )
{
	free( grid -> buffer );

	grid -> buffer = NULL;
	grid -> x = grid -> y = grid -> z = NULL;

	grid->nx = 0;
	grid->gc[0] = 0;
	grid->gc[1] = 0;

	return 0;

}

void cvfld_grid_zero( t_cvfld_grid *grid ) {
	size_t size = 3 * (grid->gc[0] + grid->nx + grid->gc[1]) * sizeof(float complex);
	memset( grid -> buffer, 0, size );
}

