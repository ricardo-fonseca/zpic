/*
 *  timer.c
 *  zpic
 *
 *  Created by Ricardo Fonseca on 13/8/10.
 *  Copyright 2010 Centro de Física dos Plasmas. All rights reserved.
 *
 */

#include "timer.h"
#include <stdlib.h>
#include <sys/time.h>

uint64_t timer_ticks()
{
	struct timeval tv;
	gettimeofday(&tv, NULL);
	
	return ((uint64_t)tv.tv_sec)*1000000 + (uint64_t)tv.tv_usec;
}

double timer_interval_seconds(uint64_t start, uint64_t end)
{
	return (end - start) * 1.0e-6;
}

double timer_cpu_seconds()
{
    struct timeval tv;
    double wtime;
    gettimeofday(&tv, NULL);
    wtime = tv.tv_sec;
    wtime += (double)tv.tv_usec * 1.0e-6;
    return wtime;
}

double timer_resolution()
{
	struct timeval tv1, tv2;
	
	gettimeofday(&tv1, NULL);
	do {
		gettimeofday(&tv2, NULL);
	} while (tv1.tv_usec == tv2.tv_usec);
	
	return (double)(tv2.tv_usec - tv1.tv_usec) * 1.0e-6;
}
