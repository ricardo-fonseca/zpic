/**
 * @file timer.c
 * @author Ricardo Fonseca
 * @brief Timing utilities
 * @version 0.1
 * @date 2022-02-04
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include "timer.h"
#include <stdlib.h>
#include <sys/time.h>

/**
 * @brief Gets current number of timer ticks
 * 
 * This implementation is based on `gettimeofday()`, so ticks correspond
 * to microsseconds since midnight (0 hour), January 1, 1970
 * 
 * @return Number of timer ticks
 */
uint64_t timer_ticks( void )
{
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return ((uint64_t)tv.tv_sec)*1000000 + (uint64_t)tv.tv_usec;
}

/**
 * @brief Time interval between to events in seconds
 * 
 * @param start 	Start event number of ticks
 * @param end 		End event number of ticks
 * @return 			Time interval in seconds between end and start
 */
double timer_interval_seconds(uint64_t start, uint64_t end)
{
	return (end - start) * 1.0e-6;
}

/**
 * @brief Gets current system seconds
 * 
 * This implementation is based on `gettimeofday()`, so it corresponds
 * to number of seconds since midnight (0 hour), January 1, 1970
 * 
 * @return double 
 */
double timer_cpu_seconds( void )
{
    struct timeval tv;
    double wtime;
    gettimeofday(&tv, NULL);
    wtime = tv.tv_sec;
    wtime += (double)tv.tv_usec * 1.0e-6;
    return wtime;
}

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
double timer_resolution()
{
	struct timeval tv1, tv2;
	gettimeofday(&tv1, NULL);
	do {
		gettimeofday(&tv2, NULL);
	} while (tv1.tv_usec == tv2.tv_usec);

	// Protection against going over to next second
	if (tv1.tv_sec < tv2.tv_sec ) tv2.tv_sec += 1000000;

	return (double)(tv2.tv_usec - tv1.tv_usec) * 1.0e-6;
}
