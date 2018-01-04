/*
 *  random.c
 *  zpic
 *
 *  Created by Ricardo Fonseca on 11/8/10.
 *  Copyright 2010 Centro de Física dos Plasmas. All rights reserved.
 *
 */

#include "random.h"
#include <math.h>


uint32_t m_w = 12345;    /* must not be zero */
uint32_t m_z = 67890;    /* must not be zero */

int iset = 1;

void set_rand_seed( uint32_t m_w_, uint32_t m_z_ )
{
	m_w = m_w_;
	m_z = m_z_;
}

uint32_t rand_uint32()
{
    m_z = 36969 * (m_z & 65535) + (m_z >> 16);
    m_w = 18000 * (m_w & 65535) + (m_w >> 16);
    return (m_z << 16) + m_w;  /* 32-bit result */
}

double rand_norm()
{
	static int iset = 0;
	static double gset = 0.0;
	
	if (iset) {
		iset = 0;
		return gset;
	} else {
		double v1, v2, rsq, fac;
		do {
			// pick 2 uniform numbers inside the square extending from -1 to 1
			// in each direction

			v1 = ( rand_uint32() + 0.5 ) / 2147483649.0 - 1.0;
			v2 = ( rand_uint32() + 0.5 ) / 2147483649.0 - 1.0;
			
			// check if they are inside the unit circle, and not (0,0) 
			// otherwise try again
			rsq = v1*v1 + v2*v2;
			
		} while ( rsq == 0.0 || rsq >= 1.0 );
		
		// Use Box-Muller method to generate random deviates with 
		// normal (gaussian) distribution
		fac = sqrt(-2.0 * log(rsq)/rsq);
		
		// store 1 value for future use
		gset = v1*fac;
		iset = 1;

		return v2*fac;
	}

}
