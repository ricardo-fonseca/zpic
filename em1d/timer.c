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

#if defined(_MSC_VER) || defined(_WIN32) || defined(_WIN64)
	// ----------------------------------------------------------------------------------------
	//	This is a drop in a replacement for the POSIX timing functions...
	// ----------------------------------------------------------------------------------------

	// a simple forward declaration.. needed in Windows for various reasons
	#include <windows.h>
		
	#if defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
		#define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64
	#else
		#define DELTA_EPOCH_IN_MICROSECS  11644473600000000ULL
	#endif

	struct timezone
	{
	  int  tz_minuteswest; /* minutes W of Greenwich */
	  int  tz_dsttime;     /* type of dst correction */
	};
	
	
	int gettimeofday(struct timeval *tv, struct timezone *tz) {
	  
	  // Define a structure to receive the current Windows filetime
	  FILETIME ft;
	 
		// Initialize the present time to 0 and the timezone to UTC
	  unsigned __int64 tmpres = 0;
	  static int tzflag = 0;
	 
	  if (NULL != tv) {
		GetSystemTimeAsFileTime(&ft);
	 
		// The GetSystemTimeAsFileTime returns the number of 100 nanosecond 
		// intervals since Jan 1, 1601 in a structure. Copy the high bits to 
		// the 64 bit tmpres, shift it left by 32 then or in the low 32 bits.
		tmpres |= ft.dwHighDateTime;
		tmpres <<= 32;
		tmpres |= ft.dwLowDateTime;
	 
		// Convert to microseconds by dividing by 10
		tmpres /= 10;
	 
		// The Unix epoch starts on Jan 1 1970.  Need to subtract the difference 
		// in seconds from Jan 1 1601.
		tmpres -= DELTA_EPOCH_IN_MICROSECS;
	 
		// Finally change microseconds to seconds and place in the seconds value. 
		// The modulus picks up the microseconds.
		tv->tv_sec = (long)(tmpres / 1000000UL);
		tv->tv_usec = (long)(tmpres % 1000000UL);
	  }
	 
	  if (NULL != tz) {
		if(!tzflag) {
			_tzset();
			tzflag++;
		}
		// Adjust for the timezone west of Greenwich
		//tz->tz_minuteswest = _timezone / 60;
		//tz->tz_dsttime = _daylight;
	  }
	  return 0;
	}
#else
	#include <sys/time.h>
#endif

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
