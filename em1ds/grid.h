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

/**
 * @brief Scalar float grid
 * 
 */
typedef struct ScalarGrid {

	float *s;		///< Pointer to position 0
	float *buffer;	///< Data buffer
	int nx;			///< Number of points (excludes guard cells)
	int gc[2];		///< Number of guard cells (lower/upper)

} t_scalar_grid;

/**
 * @brief Scalar complex float grid
 * 
 */
typedef struct ComplexScalarGrid {

	float complex *s;		///< Pointer to position 0
	float complex *buffer;	///< Data buffer
	int nx;					///< Number of points (excludes guard cells)
	int gc[2];				///< Number of guard cells (lower/upper)

} t_cscalar_grid;

/**
 * @brief 3 component vector float grid
 * 
 * Data elements are kept using a SoA (structure of arrays) arangement
 * 
 */
typedef struct Vector3Grid{

	float 	*x,		///< Pointer to position 0 of x-values grid
			*y,		///< Pointer to position 0 of y-values grid
			*z;		///< Pointer to position 0 of z-values grid
	float *buffer;	///< Data buffer
	int nx;			///< Numnber of grid points (excludes guard cells)
	int gc[2];		///< Number of guard cells (lower / upper)

} t_float3_grid;

/**
 * @brief 3 component vector complex float grid
 * 
 * Data elements are kept using a SoA (structure of arrays) arangement
 * 
 */
typedef struct ComplexVector3Grid{

	float complex 	*x,		///< Pointer to position 0 of x-values grid
					*y,		///< Pointer to position 0 of y-values grid
					*z;		///< Pointer to position 0 of z-values grid
	float complex *buffer;	///< Data buffer
	int nx;					///< Numnber of grid points (excludes guard cells)
	int gc[2];				///< Number of guard cells (lower / upper)

} t_cfloat3_grid;

/**
 * @brief Initialize ScalarGrid variable
 * 
 * @param grid 		Scalar grid
 * @param nx 		Number of points (excludes guard cells)
 * @param gc 		Number of guard cells (lower/upper), may be set to NULL
 * 					specifying 0 guard cells
 * @return			0 on success, -1 on error
 */
int  scalar_grid_init( t_scalar_grid *grid, const int nx, const int * gc );

/**
 * @brief Free dynamic memory from ScalarGrid variable
 * 
 * @param grid 	Scalar grid
 * @return 		0 on success (always returns 0)
 */
int  scalar_grid_cleanup( t_scalar_grid *grid );

/**
 * @brief Sets all grid values to zero
 * 
 * @param grid 	Scalar grid
 */
void scalar_grid_zero( t_scalar_grid *grid );

/**
 * @brief Copies values from one scalar grid to another
 * 
 * The destination grid needs to be allocated beforehand and must have the
 * same size (including guard cells) as the source grid. Also, the source
 * and destination buffers must not overlap, otherwise behavior is
 * undefined
 * 
 * @param dst 	Destination grid
 * @param src 	Source grid
 */
void scalar_grid_copy( t_scalar_grid *dst, t_scalar_grid *src );

/**
 * @brief Initialize ComplexScalarGrid variable
 * 
 * @param grid 		Complex scalar grid
 * @param nx 		Number of points (excludes guard cells)
 * @param gc 		Number of guard cells (lower/upper), may be set to NULL
 * 					specifying 0 guard cells
 * @return			0 on success, -1 on error
 */
int  cscalar_grid_init( t_cscalar_grid *grid, const int nx, const int * gc );

/**
 * @brief Free dynamic memory from ComplexScalarGrid variable
 * 
 * @param grid 	Scalar grid
 * @return 		0 on success (always returns 0)
 */
int  cscalar_grid_cleanup( t_cscalar_grid *grid );

/**
 * @brief Sets all grid values to zero
 * 
 * @param grid 	Complex scalar grid
 */
void cscalar_grid_zero( t_cscalar_grid *grid );

/**
 * @brief Initialize Vector3Grid variable
 * 
 * @param grid 		Vector3 grid
 * @param nx 		Number of points (excludes guard cells)
 * @param gc 		Number of guard cells (lower/upper), may be set to NULL
 * 					specifying 0 guard cells
 * @return			0 on success, -1 on error
 */
int  float3_grid_init( t_float3_grid *grid, const int nx, const int * gc );

/**
 * @brief Free dynamic memory from Vector3Grid variable
 * 
 * @param grid 	Vector3 grid
 * @return 		0 on success (always returns 0)
 */
int  float3_grid_cleanup( t_float3_grid *grid );

/**
 * @brief Sets all grid values to zero
 * 
 * @param grid 	Vector3 grid
 */
void float3_grid_zero( t_float3_grid *grid );

/**
 * @brief Initialize ComplexVector3Grid variable
 * 
 * @param grid 		ComplexVector3 grid
 * @param nx 		Number of points (excludes guard cells)
 * @param gc 		Number of guard cells (lower/upper), may be set to NULL
 * 					specifying 0 guard cells
 * @return			0 on success, -1 on error
 */
int  cfloat3_grid_init( t_cfloat3_grid *grid, const int nx, const int * gc );

/**
 * @brief Free dynamic memory from ComplexVector3Grid variable
 * 
 * @param grid 	ComplexVector3 grid
 * @return 		0 on success (always returns 0)
 */
int  cfloat3_grid_cleanup( t_cfloat3_grid *grid );

/**
 * @brief Sets all grid values to zero
 * 
 * @param grid 	ComplextVector3 grid
 */
void cfloat3_grid_zero( t_cfloat3_grid *grid );

#endif
