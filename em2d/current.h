/*
 *  current.h
 *  zpic
 *
 *  Created by Ricardo Fonseca on 12/8/10.
 *  Copyright 2010 Centro de Física dos Plasmas. All rights reserved.
 *
 */

#ifndef __CURRENT__
#define __CURRENT__

#include "zpic.h"

/**
 * @brief Types of digital filtering
 * 
 */
enum smooth_type { 
	NONE,		///< No filtering 
	BINOMIAL,	///< Binomial filtering
	COMPENSATED	///< Compensated binomial filtering
};

/**
 * @brief Digital filtering parameters
 * 
 * Stores digital filtering parameters
 */
typedef struct Smooth {
	enum smooth_type xtype;	///< Type of digital filtering along x
	enum smooth_type ytype;	///< Type of digital filtering along y
	int xlevel;				///< Level of filtering along x
	int ylevel;				///< Level of filtering along y
} t_smooth;

/**
 * @brief Current density object
 * 
 */
typedef struct Current {
	
	float3 *J;		///< Pointer to grid cell [0,0]
	
	float3 *J_buf;	///< Current density buffer (includes guard cells)
	
	// Grid parameters
	int nx[2];		///< Number of grid points [x,y] (excluding guard cells)
	int nrow;		///< Stride along y
	int gc[2][2];	///< Number of guard cells [x,y][lower,upper]
	
	// Box size
	float box[2];	///< Physical size of simulation box [x,y]
	
	// Cell size
	float dx[2];	///< Grid cell size [x,y]

	// Current smoothing
	t_smooth smooth;	///< Digital filtering parameters

	/// Time step
	float dt;

	/// Iteration number
	int iter;

	/// Moving window
	int moving_window;
	
} t_current;

/**
 * @brief Initializes Electric current density object
 * 
 * @param current 	Electric current density
 * @param nx 		Number of cells
 * @param box 		Physical box size
 * @param dt 		Simulation time step
  */
void current_new( t_current *current, int nx[], float box[], float dt );

/**
 * @brief Frees dynamic memory from electric current density
 * 
 * @param current Electric current density object
 */
void current_delete( t_current *current );

/**
 * @brief Sets all electric current density values to zero
 * 
 * @param current Electric current density object
 */
void current_zero( t_current *current );

/**
 * @brief Advances electric current density 1 time step
 * 
 * @param current Electric current density object
 */
void current_update( t_current *current );

/**
 * @brief Saves electric current density diagnostic information to disk
 * 
 * @param current Electric current density object
 * @param jc Current component to save, must be one of {0,1,2}
 */
void current_report( const t_current *current, const int jc );

#endif
