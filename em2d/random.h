/*
 *  random.h
 *  zpic
 *
 *  Created by Ricardo Fonseca on 11/8/10.
 *  Copyright 2010 Centro de Física dos Plasmas. All rights reserved.
 *
 */

#ifndef __RANDOM__
#define __RANDOM__

#include <stdint.h>

/**
 * @brief Sets the seed for the pseudo random number generator
 * 
 * @param m_w_ Seed value w (must not be zero)
 * @param m_z_ Seed value z (must not be zero)
 */
void set_rand_seed( uint32_t m_z_, uint32_t m_w_ );

/**
 * @brief Returns a variate of the normal distribution (mean 0, stdev 1)
 * 
 * @return Double precision random number following a normal distribution 
 */
double rand_norm( void );

#endif
