/**
 * @file grid2d.c
 * @author Ricardo Fonseca
 * @brief 2D grids
 * @version 0.2
 * @date 2022-02-04
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include "grid2d.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/****************************************************************************************
	Scalar grids
 ****************************************************************************************/

/**
 * @brief Initialize ScalarGrid2D variable
 * 
 * @param grid 		2D scalar grid
 * @param nx 		Number of points [x,y] (excludes guard cells)
 * @param gc 		Number of guard cells [x,y][lower,upper], may be set to NULL
 * @return 			0 on success, -1 on error
 */
int scalar_grid2d_init( t_scalar_grid2d *grid, const unsigned int nx[], const unsigned int gc[][2] )
{

	// store nx and gc values
	grid->nx[0] = nx[0];
	grid->nx[1] = nx[1];

	if ( gc ) {
		grid->gc[0][0] = gc[0][0];
		grid->gc[0][1] = gc[0][1];
		grid->gc[1][0] = gc[1][0];
		grid->gc[1][1] = gc[1][1];
	} else {
		grid->gc[0][0] = grid->gc[0][1] = 0;
		grid->gc[1][0] = grid->gc[1][1] = 0;
	}

	size_t size = (grid->gc[0][0] + grid->nx[0] + grid->gc[0][1]) *
				  (grid->gc[1][0] + grid->nx[1] + grid->gc[1][1]);
	
	grid -> buffer = (float *) malloc( size * sizeof( float ) );
	
	if ( !grid -> buffer ) {
		fprintf(stderr, "(*error*) Unable to allocate memory for fld_grid variable\n");
		return(-1);
	}

	grid -> nrow = grid->gc[0][0] + grid->nx[0] + grid->gc[0][1];

	// Make s point to cell [0]
	grid -> s = grid -> buffer + grid->gc[0][0] + grid->gc[1][0] * grid -> nrow ;

	return 0;
}

/**
 * @brief Free dynamic memory from ScalarGrid2D variable
 * 
 * @param grid 	2D scalar grid
 * @return 		0 on success (always returns 0)
 */
int scalar_grid2d_cleanup( t_scalar_grid2d *grid )
{
	free( grid -> buffer );

	grid -> buffer = NULL;
	grid -> s = NULL;

	grid->nx[0] = grid->nx[1] = 0;
	grid->gc[0][0] = grid->gc[0][1] = 0;
	grid->gc[1][0] = grid->gc[1][1] = 0;

	return 0;
}

/**
 * @brief Sets all grid values to zero
 * 
 * @param grid 	2D scalar grid
 */
void scalar_grid2d_zero( t_scalar_grid2d *grid ) {
	size_t size = (grid->gc[0][0] + grid->nx[0] + grid->gc[0][1]) *
				  (grid->gc[1][0] + grid->nx[1] + grid->gc[1][1]) *
				  sizeof(float);
	memset( grid -> buffer, 0, size );
}

/**
 * @brief Copies values from one 2D scalar grid to another
 * 
 * The destination grid needs to be allocated beforehand and must have the
 * same size (including guard cells) as the source grid. Also, the source
 * and destination buffers must not overlap, otherwise behavior is
 * undefined
 * 
 * @param dst 	Destination grid
 * @param src 	Source grid
 */
void scalar_grid2d_copy( t_scalar_grid2d * restrict dst, const t_scalar_grid2d * restrict src ){
	size_t size = (src->gc[0][0] + src->nx[0] + src->gc[0][1]) *
				  (src->gc[1][0] + src->nx[1] + src->gc[1][1]) *
				  sizeof(float);
	memcpy( dst->buffer, src->buffer, size);
}

/**
 * @brief Adds 2 2D scalar grids in place
 * 
 * The routine performs dst = dst + src
 * 
 * @param dst 	Destination grid
 * @param src 	Source grid
 */
void scalar_grid2d_add( t_scalar_grid2d *dst, const t_scalar_grid2d *src ){
	unsigned long size = (src->gc[0][0] + src->nx[0] + src->gc[0][1]) *
				         (src->gc[1][0] + src->nx[1] + src->gc[1][1]);

	unsigned long i;
	for( i = 0; i < size; i++) dst->buffer[i] += src->buffer[i];
}


/****************************************************************************************
	Complex Scalar grids
 ****************************************************************************************/

/**
 * @brief Initialize ComplexScalarGrid2D variable
 * 
 * @param grid 		2D complex scalar grid
 * @param nx 		Number of points [x,y] (excludes guard cells)
 * @param gc 		Number of guard cells [x,y][lower,upper], may be set to NULL
 * @return 			0 on success, -1 on error
 */
int cscalar_grid2d_init( t_cscalar_grid2d *grid, const unsigned int nx[], const unsigned int gc[][2] )
{

	// store nx and gc values
	grid->nx[0] = nx[0];
	grid->nx[1] = nx[1];

	if ( gc ) {
		grid->gc[0][0] = gc[0][0];
		grid->gc[0][1] = gc[0][1];
		grid->gc[1][0] = gc[1][0];
		grid->gc[1][1] = gc[1][1];
	} else {
		grid->gc[0][0] = grid->gc[0][1] = 0;
		grid->gc[1][0] = grid->gc[1][1] = 0;
	}


	size_t size = (grid->gc[0][0] + grid->nx[0] + grid->gc[0][1]) *
				  (grid->gc[1][0] + grid->nx[1] + grid->gc[1][1]);

	
	grid -> buffer = (float complex *) malloc( size * sizeof( float complex ) );
	
	if ( !grid -> buffer ) {
		fprintf(stderr, "(*error*) Unable to allocate memory for fld_grid variable\n");
		return(-1);
	}

	grid -> nrow = grid->gc[0][0] + grid->nx[0] + grid->gc[0][1];

	// Make s point to cell [0]
	grid -> s = grid -> buffer + grid->gc[0][0] + grid->gc[1][0] * grid -> nrow ;	

	return 0;
}

/**
 * @brief Free dynamic memory from ComplexScalarGrid2D variable
 * 
 * @param grid 	2D complex scalar grid
 * @return 		0 on success (always returns 0)
 */
int cscalar_grid2d_cleanup( t_cscalar_grid2d *grid )
{
	free( grid -> buffer );

	grid -> buffer = NULL;
	grid -> s = NULL;

	grid->nx[0] = grid->nx[1] = 0;
	grid->gc[0][0] = grid->gc[0][1] = 0;
	grid->gc[1][0] = grid->gc[1][1] = 0;

	return 0;
}

/**
 * @brief Sets all grid values to zero
 * 
 * @param grid 	2D complex scalar grid
 */
void cscalar_grid2d_zero( t_cscalar_grid2d *grid ) {
	size_t size = (grid->gc[0][0] + grid->nx[0] + grid->gc[0][1]) *
				  (grid->gc[1][0] + grid->nx[1] + grid->gc[1][1]) * 
				  sizeof(float complex);
	memset( grid -> buffer, 0, size );
}

/**
 * @brief Copies values from one 2D complex scalar grid to another
 * 
 * The destination grid needs to be allocated beforehand and must have the
 * same size (including guard cells) as the source grid. Also, the source
 * and destination buffers must not overlap, otherwise behavior is
 * undefined
 * 
 * @param dst 	Destination grid
 * @param src 	Source grid
 */
void cscalar_grid2d_copy( t_cscalar_grid2d * restrict dst, const t_cscalar_grid2d * restrict src ){
	size_t size = (src->gc[0][0] + src->nx[0] + src->gc[0][1]) *
				  (src->gc[1][0] + src->nx[1] + src->gc[1][1]) *
				  sizeof(float complex);
	memcpy( dst->buffer, src->buffer, size);
}

/**
 * @brief Adds 2 2D complex scalar grids in place
 * 
 * The routine performs dst = dst + src
 * 
 * @param dst 	Destination grid
 * @param src 	Source grid
 */
void cscalar_grid2d_add( t_cscalar_grid2d *dst, const t_cscalar_grid2d *src ){
	unsigned long size = (src->gc[0][0] + src->nx[0] + src->gc[0][1]) *
				         (src->gc[1][0] + src->nx[1] + src->gc[1][1]);

	unsigned long i;
	for( i = 0; i < size; i++) dst->buffer[i] += src->buffer[i];
}

/****************************************************************************************
	Vector field grids
 ****************************************************************************************/

/**
 * @brief Initialize Vector3Grid2D variable
 * 
 * @param grid 		2D vector3 grid
 * @param nx 		Number of points [x,y] (excludes guard cells)
 * @param gc 		Number of guard cells [x,y][lower,upper], may be set to NULL
 * @return 			0 on success, -1 on error
 */
int float3_grid2d_init( t_float3_grid2d *grid, const unsigned int nx[], const unsigned int gc[][2] )
{

	// store nx and gc values
	grid->nx[0] = nx[0];
	grid->nx[1] = nx[1];

	if ( gc ) {
		grid->gc[0][0] = gc[0][0];
		grid->gc[0][1] = gc[0][1];
		grid->gc[1][0] = gc[1][0];
		grid->gc[1][1] = gc[1][1];
	} else {
		grid->gc[0][0] = grid->gc[0][1] = 0;
		grid->gc[1][0] = grid->gc[1][1] = 0;
	}

	size_t size = (grid->gc[0][0] + grid->nx[0] + grid->gc[0][1]) *
				  (grid->gc[1][0] + grid->nx[1] + grid->gc[1][1]);
	
	grid -> buffer = malloc( 3 * size * sizeof( float ) );
	
	if ( !grid -> buffer ) {
		fprintf(stderr, "(*error*) Unable to allocate memory for fld_grid variable\n");
		return(-1);
	}

	grid -> nrow = grid->gc[0][0] + grid->nx[0] + grid->gc[0][1];
	
	// Make x, y and z point to cell [0]
	grid -> x = grid -> buffer + grid->gc[0][0] + grid->gc[1][0] * grid -> nrow;
	grid -> y = grid -> x + size;
	grid -> z = grid -> y + size;

	return 0;
}

/**
 * @brief Free dynamic memory from Vector3Grid2D variable
 * 
 * @param grid 	2D vector3 grid
 * @return 		0 on success (always returns 0)
 */
int float3_grid2d_cleanup( t_float3_grid2d *grid )
{
	free( grid -> buffer );

	grid -> buffer = NULL;
	grid -> x = grid -> y = grid -> z = NULL;
	grid->nx[0] = grid->nx[1] = 0;
	grid->gc[0][0] = grid->gc[0][1] = 0;
	grid->gc[1][0] = grid->gc[1][1] = 0;

	return 0;
}

/**
 * @brief Sets all grid values to zero
 * 
 * @param grid 	2D float3 grid
 */
void float3_grid2d_zero( t_float3_grid2d *grid ) {
	size_t size = (grid->gc[0][0] + grid->nx[0] + grid->gc[0][1]) *
				  (grid->gc[1][0] + grid->nx[1] + grid->gc[1][1]) *
				  3 * sizeof(float);
	memset( grid -> buffer, 0, size );
}

/**
 * @brief Copies values from one 2D scalar grid to another
 * 
 * The destination grid needs to be allocated beforehand and must have the
 * same size (including guard cells) as the source grid. Also, the source
 * and destination buffers must not overlap, otherwise behavior is
 * undefined
 * 
 * @param dst 	Destination grid
 * @param src 	Source grid
 */
void float3_grid2d_copy( t_float3_grid2d * restrict dst, const t_float3_grid2d * restrict src ){
	size_t size = (src->gc[0][0] + src->nx[0] + src->gc[0][1]) *
				  (src->gc[1][0] + src->nx[1] + src->gc[1][1]) *
				  3 * sizeof(float);
	memcpy( dst->buffer, src->buffer, size);
}

/**
 * @brief Adds 2 2D float3 grids in place
 * 
 * The routine performs dst = dst + src
 * 
 * @param dst 	Destination grid
 * @param src 	Source grid
 */
void float3_grid2d_add( t_float3_grid2d * restrict dst, const t_float3_grid2d * restrict src ){
	unsigned long size = (src->gc[0][0] + src->nx[0] + src->gc[0][1]) *
				         (src->gc[1][0] + src->nx[1] + src->gc[1][1]) * 3;

	unsigned long i;
	for( i = 0; i < size; i++) dst->buffer[i] += src->buffer[i];
}

/****************************************************************************************
	Complex Vector field grids
 ****************************************************************************************/

/**
 * @brief Initialize ComplexVector3Grid2D variable
 * 
 * @param grid 		2D complex vector3 grid
 * @param nx 		Number of points [x,y] (excludes guard cells)
 * @param gc 		Number of guard cells [x,y][lower,upper], may be set to NULL
 * @return 			0 on success, -1 on error
 */
int cfloat3_grid2d_init( t_cfloat3_grid2d *grid, const unsigned int nx[], const unsigned int gc[][2] )
{

	// store nx and gc values
	grid->nx[0] = nx[0];
	grid->nx[1] = nx[1];

	if ( gc ) {
		grid->gc[0][0] = gc[0][0];
		grid->gc[0][1] = gc[0][1];
		grid->gc[1][0] = gc[1][0];
		grid->gc[1][1] = gc[1][1];
	} else {
		grid->gc[0][0] = grid->gc[0][1] = 0;
		grid->gc[1][0] = grid->gc[1][1] = 0;
	}

	size_t size = (grid->gc[0][0] + grid->nx[0] + grid->gc[0][1]) *
				  (grid->gc[1][0] + grid->nx[1] + grid->gc[1][1]);
	
	grid -> buffer = malloc( 3 * size * sizeof( float complex ) );
	
	if ( !grid -> buffer ) {
		fprintf(stderr, "(*error*) Unable to allocate memory for fld_grid variable\n");
		return(-1);
	}

	grid -> nrow = grid->gc[0][0] + grid->nx[0] + grid->gc[0][1];
	
	// Make x, y and z point to cell [0]
	grid -> x = grid -> buffer + grid->gc[0][0] + grid->gc[1][0] * grid -> nrow;
	grid -> y = grid -> x + size;
	grid -> z = grid -> y + size;

	return 0;
}

/**
 * @brief Free dynamic memory from ComplexVector3Grid2D variable
 * 
 * @param grid 	2D complex vector3 grid
 * @return 		0 on success (always returns 0)
 */
int cfloat3_grid2d_cleanup( t_cfloat3_grid2d *grid )
{
	free( grid -> buffer );

	grid -> buffer = NULL;
	grid -> x = grid -> y = grid -> z = NULL;

	grid->nx[0] = grid->nx[1] = 0;
	grid->gc[0][0] = grid->gc[0][1] = 0;
	grid->gc[1][0] = grid->gc[1][1] = 0;

	return 0;

}

/**
 * @brief Sets all grid values to zero
 * 
 * @param grid 	2D complex float3 grid
 */
void cfloat3_grid2d_zero( t_cfloat3_grid2d *grid ) {
	size_t size =  (grid->gc[0][0] + grid->nx[0] + grid->gc[0][1]) *
				   (grid->gc[1][0] + grid->nx[1] + grid->gc[1][1]) *
				   3 * sizeof(float complex);
	memset( grid -> buffer, 0, size );
}

/**
 * @brief Copies values from one 2D scalar grid to another
 * 
 * The destination grid needs to be allocated beforehand and must have the
 * same size (including guard cells) as the source grid. Also, the source
 * and destination buffers must not overlap, otherwise behavior is
 * undefined
 * 
 * @param dst 	Destination grid
 * @param src 	Source grid
 */
void cfloat3_grid2d_copy( t_cfloat3_grid2d * restrict dst, const t_cfloat3_grid2d * restrict src ){
	size_t size = (src->gc[0][0] + src->nx[0] + src->gc[0][1]) *
				  (src->gc[1][0] + src->nx[1] + src->gc[1][1]) *
				  3 * sizeof(float complex);
	memcpy( dst->buffer, src->buffer, size);
}

/**
 * @brief Adds 2 2D complex float3 grids in place
 * 
 * The routine performs dst = dst + src
 * 
 * @param dst 	Destination grid
 * @param src 	Source grid
 */
void cfloat3_grid2d_add( t_cfloat3_grid2d * restrict dst, const t_cfloat3_grid2d * restrict src ){
	unsigned long size = (src->gc[0][0] + src->nx[0] + src->gc[0][1]) *
				         (src->gc[1][0] + src->nx[1] + src->gc[1][1]) * 3;

	unsigned long i;
	for( i = 0; i < size; i++) dst->buffer[i] += src->buffer[i];
}


