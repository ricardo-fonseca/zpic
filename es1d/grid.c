/**
 * @file grid.c
 * @author Ricardo Fonseca
 * @brief Scalar grids
 * @version 0.2
 * @date 2022-02-21
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include "grid.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/****************************************************************************************
	Scalar grids
 ****************************************************************************************/

/**
 * @brief Initialize ScalarGrid variable
 * 
 * @param grid 		Scalar grid
 * @param nx 		Number of points (excludes guard cells)
 * @param gc 		Number of guard cells (lower/upper), may be set to NULL
 * 					specifying 0 guard cells
 * @return			0 on success, -1 on error
 */
int scalar_grid_init( t_scalar_grid *grid, const int nx, const int * gc )
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

/**
 * @brief Free dynamic memory from ScalarGrid variable
 * 
 * @param grid 	Scalar grid
 * @return 		0 on success (always returns 0)
 */
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

/**
 * @brief Sets all grid values to zero
 * 
 * @param grid 	Scalar grid
 */
void scalar_grid_zero( t_scalar_grid *grid ) {
	size_t size =  (grid->gc[0] + grid->nx + grid->gc[1]) * sizeof(float);
	memset( grid -> buffer, 0, size );
}

/**
 * @brief Copies values from one scalar grid to another
 * 
 * The destination grid needs to be allocated beforehand and must have the
 * same size (including guard cells) as the source grid. Also the source
 * and destination buffers must not overlap, otherwise behavior is
 * undefined
 * 
 * @param dst 	Destination grid
 * @param src 	Source grid
 */
void scalar_grid_copy( t_scalar_grid *dst, t_scalar_grid *src ){
	size_t n = ( dst -> gc[0] + dst -> nx + dst -> gc[1] ) * sizeof( float );
	memcpy( dst->buffer, src->buffer, n);
}


/****************************************************************************************
	Complex Scalar grids
 ****************************************************************************************/

/**
 * @brief Initialize ComplexScalarGrid variable
 * 
 * @param grid 		Complex scalar grid
 * @param nx 		Number of points (excludes guard cells)
 * @param gc 		Number of guard cells (lower/upper), may be set to NULL
 * 					specifying 0 guard cells
 * @return			0 on success, -1 on error
 */
int cscalar_grid_init( t_cscalar_grid *grid, const int nx, const int * gc )
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

/**
 * @brief Free dynamic memory from ComplexScalarGrid variable
 * 
 * @param grid 	Scalar grid
 * @return 		0 on success (always returns 0)
 */
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

/**
 * @brief Sets all grid values to zero
 * 
 * @param grid 	Complex scalar grid
 */
void cscalar_grid_zero( t_cscalar_grid *grid ) {
	size_t size =  (grid->gc[0] + grid->nx + grid->gc[1]) * sizeof(float complex);
	memset( grid -> buffer, 0, size );
}
