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

/**
 * @brief Scalar float 2D grid
 * 
 */
typedef struct ScalarGrid2D {

	float *s;		///< Pointer to position [0,0]

	float *buffer;	///< Data buffer

	int nx[2];		///< Number of points [x,y] (excludes guard cells)
	int nrow;		///< Stride along y
	int gc[2][2];	///< Number of guard cells [x,y][lower,upper]

} t_scalar_grid2d;

/**
 * @brief Scalar complex float 2D grid
 * 
 */
typedef struct ComplexScalarGrid2D {

	float complex *s;		///< Pointer to position [0,0]

	float complex *buffer;	///< Data buffer

	int nx[2];		///< Number of points [x,y] (excludes guard cells)
	int nrow;		///< Stride along y
	int gc[2][2];	///< Number of guard cells [x,y][lower,upper]

} t_cscalar_grid2d;

/**
 * @brief 3 component vector float 2D grid
 * 
 * Data elements are kept using a SoA (structure of arrays) arangement
 * 
 */
typedef struct Vector3Grid2D {

	float 	*x,	///< Pointer to position [0,0] of x-values grid
			*y,	///< Pointer to position [0,0] of y-values grid
			*z;	///< Pointer to position [0,0] of z-values grid

	float *buffer; ///< Data buffer

	int nx[2];		///< Number of points [x,y] (excludes guard cells)
	int nrow;		///< Stride along y
	int gc[2][2];	///< Number of guard cells [x,y][lower,upper]
} t_float3_grid2d;

/**
 * @brief 3 component vector complex float 2D grid
 * 
 * Data elements are kept using a SoA (structure of arrays) arangement
 * 
 */
typedef struct ComplexVector3Grid2D {

	float complex 	*x,	///< Pointer to position [0,0] of x-values grid
					*y,	///< Pointer to position [0,0] of y-values grid
					*z;	///< Pointer to position [0,0] of z-values grid

	float complex *buffer;	///< Data buffer

	int nx[2];		///< Number of points [x,y] (excludes guard cells)
	int nrow;		///< Stride along y
	int gc[2][2];	///< Number of guard cells [x,y][lower,upper]

} t_cfloat3_grid2d;

/**
 * @brief Initialize ScalarGrid2D variable
 * 
 * @param grid 		2D scalar grid
 * @param nx 		Number of points [x,y] (excludes guard cells)
 * @param gc 		Number of guard cells [x,y][lower,upper], may be set to NULL
 * @return 			0 on success, -1 on error
 */
int  scalar_grid2d_init( t_scalar_grid2d *grid, const unsigned int nx[], const unsigned int gc[][2] );

/**
 * @brief Free dynamic memory from ScalarGrid2D variable
 * 
 * @param grid 	2D scalar grid
 * @return 		0 on success (always returns 0)
 */
int  scalar_grid2d_cleanup( t_scalar_grid2d *grid );

/**
 * @brief Sets all grid values to zero
 * 
 * @param grid 	2D scalar grid
 */
void scalar_grid2d_zero( t_scalar_grid2d *grid );

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
void scalar_grid2d_copy( t_scalar_grid2d * restrict dst, const t_scalar_grid2d * restrict src );

/**
 * @brief Adds 2 2D scalar grids in place
 * 
 * The routine performs dst = dst + src
 * 
 * @param dst 	Destination grid
 * @param src 	Source grid
 */
void scalar_grid2d_add( t_scalar_grid2d * restrict dst, const t_scalar_grid2d * restrict src );

/**
 * @brief Initialize ComplexScalarGrid2D variable
 * 
 * @param grid 		2D complex scalar grid
 * @param nx 		Number of points [x,y] (excludes guard cells)
 * @param gc 		Number of guard cells [x,y][lower,upper], may be set to NULL
 * @return 			0 on success, -1 on error
 */
int  cscalar_grid2d_init( t_cscalar_grid2d *grid, const unsigned int nx[], const unsigned int gc[][2] );

/**
 * @brief Free dynamic memory from ComplexScalarGrid2D variable
 * 
 * @param grid 	2D complex scalar grid
 * @return 		0 on success (always returns 0)
 */
int  cscalar_grid2d_cleanup( t_cscalar_grid2d *grid );

/**
 * @brief Sets all grid values to zero
 * 
 * @param grid 	2D complex scalar grid
 */
void cscalar_grid2d_zero( t_cscalar_grid2d *grid );

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
void cscalar_grid2d_copy( t_cscalar_grid2d * restrict dst, const t_cscalar_grid2d * restrict src );

/**
 * @brief Adds 2 2D complex scalar grids in place
 * 
 * The routine performs dst = dst + src
 * 
 * @param dst 	Destination grid
 * @param src 	Source grid
 */
void cscalar_grid2d_add( t_cscalar_grid2d * restrict dst, const t_cscalar_grid2d * restrict src );

/**
 * @brief Initialize Vector3Grid2D variable
 * 
 * @param grid 		2D vector3 grid
 * @param nx 		Number of points [x,y] (excludes guard cells)
 * @param gc 		Number of guard cells [x,y][lower,upper], may be set to NULL
 * @return 			0 on success, -1 on error
 */
int  float3_grid2d_init( t_float3_grid2d *grid, const unsigned int nx[], const unsigned int gc[][2] );

/**
 * @brief Free dynamic memory from Vector3Grid2D variable
 * 
 * @param grid 	2D vector3 grid
 * @return 		0 on success (always returns 0)
 */
int  float3_grid2d_cleanup( t_float3_grid2d *grid );

/**
 * @brief Sets all grid values to zero
 * 
 * @param grid 	2D float3 grid
 */
void float3_grid2d_zero( t_float3_grid2d *grid );

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
void float3_grid2d_copy( t_float3_grid2d * restrict dst, const t_float3_grid2d * restrict src );

/**
 * @brief Adds 2 2D float3 grids in place
 * 
 * The routine performs dst = dst + src
 * 
 * @param dst 	Destination grid
 * @param src 	Source grid
 */
void float3_grid2d_add( t_float3_grid2d * restrict dst, const t_float3_grid2d * restrict src );

/**
 * @brief Initialize ComplexVector3Grid2D variable
 * 
 * @param grid 		2D complex vector3 grid
 * @param nx 		Number of points [x,y] (excludes guard cells)
 * @param gc 		Number of guard cells [x,y][lower,upper], may be set to NULL
 * @return 			0 on success, -1 on error
 */
int  cfloat3_grid2d_init( t_cfloat3_grid2d *grid, const unsigned int nx[], const unsigned int gc[][2] );

/**
 * @brief Free dynamic memory from ComplexVector3Grid2D variable
 * 
 * @param grid 	2D complex vector3 grid
 * @return 		0 on success (always returns 0)
 */
int  cfloat3_grid2d_cleanup( t_cfloat3_grid2d *grid );

/**
 * @brief Sets all grid values to zero
 * 
 * @param grid 	2D complex float3 grid
 */
void cfloat3_grid2d_zero( t_cfloat3_grid2d *grid );

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
void cfloat3_grid2d_copy( t_cfloat3_grid2d * restrict dst, const t_cfloat3_grid2d * restrict src );

/**
 * @brief Adds 2 2D complex float3 grids in place
 * 
 * The routine performs dst = dst + src
 * 
 * @param dst 	Destination grid
 * @param src 	Source grid
 */
void cfloat3_grid2d_add( t_cfloat3_grid2d * restrict dst, const t_cfloat3_grid2d * restrict src );

#endif
