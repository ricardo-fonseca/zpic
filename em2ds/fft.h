/*
 *  fft.h
 *
 */

#ifndef __FFT__
#define __FFT__

/*
We should use the "Type generic math <tgmath.h> instead of <math.h> and <complex.h>
And then change the code so it can be easily compiled in single our double precision.
*/


#include <complex.h>

#ifndef M_PI
#define M_PI        3.14159265358979323846264338327950288   /* pi             */
#endif

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


typedef struct {
	t_fft_cfg cfg;
    float complex *phase;
} t_fftr_cfg;

typedef struct {
    unsigned int nx;
    unsigned int ny;
		unsigned int stridey;
    t_fftr_cfg cfgx;
    t_fft_cfg  cfgy;
} t_fftr2d_cfg;

/**********************************************************************************
Complex data FFT routines
***********************************************************************************/

void fft_c2c( t_fft_cfg* cfg, const float complex in[], float complex out[] );
int fft_init_cfg( t_fft_cfg* cfg, unsigned int n, enum fft_direction direction );
int fft_cleanup_cfg( t_fft_cfg* cfg );
float fft_dk( const unsigned int n, const float dx);


/**********************************************************************************
Real data FFT routines
***********************************************************************************/

int fftr_init_cfg( t_fftr_cfg* rcfg, unsigned int nr, enum fft_direction direction );
int fftr_cleanup_cfg( t_fftr_cfg* rcfg );

void fftr_r2c( t_fftr_cfg* rcfg, const float in[], float complex out[] );
void fftr_c2r( t_fftr_cfg* rcfg, const float complex in[], float out[] );

/**********************************************************************************
2D Real data FFT routines
***********************************************************************************/

int fftr2d_init_cfg( t_fftr2d_cfg* rcfg, unsigned int nx, unsigned int ny,
		unsigned int stridey, enum fft_direction direction );
int fftr2d_cleanup_cfg( t_fftr2d_cfg* rcfg );

void fftr2d_r2c(t_fftr2d_cfg * rcfg, const float* in, float complex * out );
void fftr2d_c2r(t_fftr2d_cfg * rcfg, const float complex * in, float *out );

#endif
