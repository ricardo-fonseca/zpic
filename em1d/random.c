/**
 * @file random.c
 * @author Ricardo Fonseca
 * @brief Pseudo-random number generator
 * @version 0.1
 * @date 2022-02-04
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include "random.h"
#include <math.h>


uint32_t m_w = 12345;    ///< Random seed w, must not be zero nor 0x464fffff
uint32_t m_z = 67890;    ///< Random seed z, must not be zero nor 0x9068ffff

/**
 * @brief Sets the seed for the pseudo random number generator
 * 
 * @param m_w_ Seed value w (must not be zero)
 * @param m_z_ Seed value z (must not be zero)
 */
void set_rand_seed( uint32_t m_w_, uint32_t m_z_ )
{
	m_w = m_w_;
	m_z = m_z_;
}

/**
 * @brief Returns a 32 bit pseudo random number using Marsaglia MWC algorithm
 * 
 * Follows George Marsaglia's post to Sci.Stat.Math on 1/12/99:
 * 
 * The  MWC generator concatenates two 16-bit multiply-
 * with-carry generators, x(n)=36969x(n-1)+carry,
 * y(n)=18000y(n-1)+carry  mod 2^16, has period about
 * 2^60 and seems to pass all tests of randomness. A
 * favorite stand-alone generator---faster than KISS,
 * which contains it. [sic]
 * 
 * See for example Numerical Recipes, 3rd edition, Section 7.1.7,
 * "When You Have Only 32-Bit Arithmetic"
 * 
 * @return Random value in the range [0,2^32 - 1]
 */
uint32_t rand_uint32( void )
{
    m_z = 36969 * (m_z & 65535) + (m_z >> 16);
    m_w = 18000 * (m_w & 65535) + (m_w >> 16);
    return (m_z << 16) + m_w;  /* 32-bit result */
}

/**
 * @brief Returns a variate of the normal distribution (mean 0, stdev 1)
 * 
 * Uses the Box-Muller method for generating random deviates with a normal
 * (Gaussian) distribution:
 * 
 *  Box, G. E. P.; Muller, Mervin E. (1958). "A Note on the Generation of 
 *  Random Normal Deviates". The Annals of Mathematical Statistics. 29 (2):
 *  610–611. doi:10.1214/aoms/1177706645.
 * 
 * @return Double precision random number following a normal distribution 
 */
double rand_norm( void )
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
