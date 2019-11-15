/*
 *  fft.h
 *
 */

#ifndef __FFT__
#define __FFT__

#include "zpic.h"
#include <complex.h>

#define MAX_FACTORS 32

enum fft_direction { FFT_FORWARD, FFT_BACKWARD };


typedef struct {
	unsigned int p, n;
} t_fft_factors;

typedef struct {
	unsigned int n;
    enum fft_direction direction;
    t_fft_factors factors[MAX_FACTORS];
    float complex *phase;
} t_fft_cfg;


enum fftr_storage { FFTR_CCS, FFTR_PERM };

typedef struct {
	enum fftr_storage storage;
	t_fft_cfg cfg;
    float complex *phase;
} t_fftr_cfg;

float fft_dk( const unsigned int n, const float dx);


/**********************************************************************************
Complex data FFT routines
***********************************************************************************/

void fft_c2c( t_fft_cfg* cfg, const float complex in[], float complex out[] );
int fft_init_cfg( t_fft_cfg* cfg, unsigned int n, enum fft_direction direction );
int fft_cleanup_cfg( t_fft_cfg* cfg );

/**********************************************************************************
Real data FFT routines
***********************************************************************************/

int fftr_init_cfg( t_fftr_cfg* rcfg, unsigned int nr, enum fft_direction direction,
    enum fftr_storage storage  );
int fftr_cleanup_cfg( t_fftr_cfg* rcfg );

void fftr_r2c( t_fftr_cfg* rcfg, const float in[], float complex out[] );
void fftr_c2r( t_fftr_cfg* rcfg, const float complex in[], float out[] );


#endif
