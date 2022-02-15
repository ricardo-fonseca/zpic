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

/**
 * @brief Gets current number of timer ticks
 * 
 * @return Number of timer ticks
 */
uint64_t timer_ticks( void );

/**
 * @brief Time interval between to events in seconds
 * 
 * @param start 	Start event number of ticks
 * @param end 		End event number of ticks
 * @return 			Time interval in seconds between end and start
 */
double timer_interval_seconds(uint64_t start, uint64_t end);

/**
 * @brief Gets current system seconds
 * 
 * @return double 
 */
double timer_cpu_seconds( void );

/**
 * @brief Gets timer resolution in seconds
 * 
 * This implementation is based on `gettimeofday()`, so the minimum
 * resolution value will be 1 microssecond.
 * 
 * Resolution is determined empirically; the routine will do successive
 * calls to `gettimeofday()` until the return time is different
 * 
 * @return Timer resolution 
 */
double timer_resolution( void );

#endif
