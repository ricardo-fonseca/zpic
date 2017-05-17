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

void set_rand_seed( uint32_t m_z_, uint32_t m_w_ );

double rand_norm();

#endif
