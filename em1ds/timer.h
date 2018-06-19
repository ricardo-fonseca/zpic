/*
 *  timer.h
 *  zpic
 *
 *  Created by Ricardo Fonseca on 13/8/10.
 *  Copyright 2010 Centro de Física dos Plasmas. All rights reserved.
 *
 */

#ifndef __TIMER__
#define __TIMER__


#include <stdint.h>

uint64_t timer_ticks( void );
double timer_interval_seconds(uint64_t start, uint64_t end);
double timer_cpu_seconds( void );
double timer_resolution( void );

#endif
