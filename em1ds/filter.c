/**
 * @file filter.c
 * @author Ricardo Fonseca
 * @brief Spectral filtering
 * @version 0.2
 * @date 2022-02-17
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#include "filter.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


#ifndef M_LN2
#define M_LN2       0.693147180559945309417232121458176568  ///< loge(2)
#endif

/**
 * @brief Updates filter values for the specified parameters
 * 
 * @param filter 	Filter object
 * @param type   	Type of filter to use
 * @param ck     	Cutoff wavenumber for the filter (filter response = 1/2)
 * @return			0 on success, -1 on error
 */
int filter_set( t_filter* const filter, enum filter_type const type, float const ck )
{

	switch( type ) {
		case( FILTER_NONE ):
			for (int i=0; i < filter->nk; i++) filter -> Sk[i] = 1.0;
			break;

		case( FILTER_GAUSS ):
		    if ( ck <= 0.0f ) {
		    	fprintf(stderr, "(*error*) Invalid valid for cutoff wavenumber, must be in the ]0.0, 1.0[ range\n");
		    	return -1;
		    }
			const float a = M_LN2 / (ck*ck);
			for (int i=0; i < filter->nk; i++) {
		        const float kx = ((float) i) / (filter->nk-1);
				filter -> Sk[i] = exp( -kx*kx * a );
			}
			break;

		case( FILTER_SHARP ):
		    if ( ck <= 0.0f || ck >= 1.0f ) {
		    	fprintf(stderr, "(*error*) Invalid valid for cutoff wavenumber, must be in the ]0.0, 1.0[ range\n");
		    	return -1;
		    }
			for (int i=0; i < filter->nk; i++) {
		        filter -> Sk[i] = ( i <= ck * filter -> nk ) ? 1.0 : 0.0;
			}
			break;

    	default:
	    	fprintf(stderr, "(*error*) Invalid valid filter type\n");
	    	return -1;

	}

	// Store the filter type
	filter -> type = type;

	return 0;
}


/**
 * @brief Initializes filter object
 * 
 * @param filter Filter object
 * @param nk     Number of k-space points (assumes k at nk-1 is Nyquist wavenumber)
 */
void filter_new( t_filter *filter, int nk ) {

	filter -> nk = nk;
	filter -> Sk = malloc( filter -> nk * sizeof(float) );
	filter_set( filter, FILTER_NONE, 0.0 );
}

/**
 * @brief Cleanup filter object
 * 
 * @param filter Filter object
 */
void filter_delete( t_filter *filter ) {
	filter -> type = FILTER_NONE;
	free( filter -> Sk );
	filter -> Sk = NULL;
}
