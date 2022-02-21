/**
 * @file fft.c
 * @author Ricardo Fonseca
 * @brief 1D/2D FFT library
 * @version 0.1
 * @date 2022-02-17
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include "fft.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/***************************************************************************************
 1D FFT of complex data
****************************************************************************************/

/**
 * @brief Performs a DFT of size 2
 * 
 * @param out       DFT input/output
 * @param stride    Data access stride
 * @param phase     Input data
 * @param m         Number of points
 */
void butterfly2( float complex * restrict const out, const unsigned int stride,
    const float complex * restrict const phase, const unsigned int m)
{
    for( unsigned int i = 0; i < m; i++) {
        float complex t0 = out[i    ];
        float complex t1 = out[i + m] * phase[i*stride];
        out[i    ] = t0 + t1;
        out[i + m] = t0 - t1;
    }
}

/**
 * @brief Performs a DFT of size 3
 * 
 * @param out       DFT input/output
 * @param stride    Data access stride
 * @param phase     Input data
 * @param m         Number of points
 */
void butterfly3( float complex * restrict const out, const unsigned int stride,
    const float complex * restrict const phase, const unsigned int m)
{
    const float complex ima = I * cimag(phase[ stride * m ]);

    for( unsigned int i = 0; i<m; i++) {
        float complex t0 = out[i      ];
        float complex t1 = out[i +   m] * phase[     i * stride ];
        float complex t2 = out[i + 2*m] * phase[ 2 * i * stride ];

        float complex s0 = t1 + t2;
        float complex s1 = t1 - t2;

        out[i      ] = t0 + s0;
        out[i +   m] = t0 - 0.5f * s0 + ima * s1;
        out[i + 2*m] = t0 - 0.5f * s0 - ima * s1;
    }
}

/**
 * @brief Performs a DFT of size 4
 * 
 * @param out       DFT input/output
 * @param stride    Data access stride
 * @param phase     Input data
 * @param m         Number of points 
 * @param direction Direction of transform {FFT_BACKWARD, FFT_FORWARD}
 */
void butterfly4( float complex * restrict const out, const unsigned int stride,
    const float complex * restrict const phase, const unsigned int m,
    const enum fft_direction direction )
{
    const float complex p = ( direction == FFT_BACKWARD )? I : -I;

    for( unsigned int i = 0; i < m; i++){
        float complex t0 = out[i      ];
        float complex t1 = out[i +   m] * phase[     i * stride ];
        float complex t2 = out[i + 2*m] * phase[ 2 * i * stride ];
        float complex t3 = out[i + 3*m] * phase[ 3 * i * stride ];

        float complex s0 = t1 + t3;
        float complex s1 = t1 - t3;
        float complex s2 = t0 + t2;
        float complex s3 = t0 - t2;

        out[i      ] = s2 + s0;
        out[i +   m] = s3 + p * s1;
        out[i + 2*m] = s2 - s0;
        out[i + 3*m] = s3 - p * s1;
    }
}

/**
 * @brief Performs a DFT of size 5
 * 
 * @param out       DFT input/output
 * @param stride    Data access stride
 * @param phase     Input data
 * @param m         Number of points
 */
void butterfly5( float complex * restrict const out, const unsigned int stride,
    const float complex * restrict const phase, const unsigned int m)
{
    const float ra = creal(phase[     stride * m ]);
    const float ia = cimag(phase[     stride * m ]);

    const float rb = creal(phase[ 2 * stride * m ]);
    const float ib = cimag(phase[ 2 * stride * m ]);

    for( unsigned int i = 0; i<m; i++) {

        float complex t0 = out[       i ];
        float complex t1 = out[   m + i ] * phase[     i * stride ];
        float complex t2 = out[ 2*m + i ] * phase[ 2 * i * stride ];
        float complex t3 = out[ 3*m + i ] * phase[ 3 * i * stride ];
        float complex t4 = out[ 4*m + i ] * phase[ 4 * i * stride ];

        float complex s0 = t1 + t4;
        float complex s1 = t1 - t4;
        float complex s2 = t2 + t3;
        float complex s3 = t2 - t3;

        float complex s4 = t0 + s0 * ra + s2 * rb;
        float complex s5 = t0 + s0 * rb + s2 * ra;

        float complex s6 = -I * (s1 * ia + s3 * ib);
        float complex s7 =  I * (s1 * ib - s3 * ia);

        out[ i       ] = t0 + s0 + s2;
        out[ i +   m ] = s4 - s6;
        out[ i + 2*m ] = s5 + s7;
        out[ i + 3*m ] = s5 - s7;
        out[ i + 4*m ] = s4 + s6;
    }
}

/**
 * @brief Performs a DFT of arbitrary size N
 * 
 * @param out       DFT input/output
 * @param stride    Data access stride
 * @param phase     Input data
 * @param m         Number of points
 * @param p         
 * @param n         
 */
void butterflyN(float complex * restrict const out, const unsigned int stride,
    const float complex * restrict const phase, const unsigned int m, const unsigned int p,
    const unsigned int n)
{
    float complex t[p];

    for( unsigned int i = 0; i<m; i++) {

        for( unsigned int j = 0; j < p; j++) {
            t[j] = out[i + m*j];
        }

        for( unsigned int j=0; j <p; j++ ) {
            unsigned int tstride = (i + j*m) * stride;
            float complex s = t[0];
            for( unsigned int k = 1; k<p; k++){
                s += t[k] * phase[ (k*tstride) % n ];
            }
            out[i + j*m] = s;
        }
    }

}

/**
 * @brief Recursively perform a Complex-to-complex FFT
 * 
 * @param in            Input data
 * @param in_stride     Input data stride
 * @param out           Output data
 * @param out_stride    Output data stride
 * @param factors       Factor decomposition of input data size
 * @param cfg           Transform configuration
 */
void transform( const float complex in[], const unsigned int in_stride,
                float complex out[],  const unsigned int out_stride,
                t_fft_factors factors[], t_fft_cfg *cfg )

{
    int p = factors[0].p;
    int m = factors[0].n;

    if ( m == 1 ) {
        for( int i = 0; i<p; i++) {
            out[i] = in[ i * in_stride * out_stride ];
        }
    } else {
        for( int i = 0; i<p; i++) {
            transform( &in[i * out_stride * in_stride], in_stride,
                       &out[i*m], out_stride * p,
                       &factors[1], cfg);
        }
    }

    switch (p){
        case 2 :
            butterfly2( out, out_stride, cfg -> phase, m ); break;
        case 3 :
            butterfly3( out, out_stride, cfg -> phase, m ); break;
        case 4 :
            butterfly4( out, out_stride, cfg -> phase, m, cfg -> direction ); break;
        case 5 :
            butterfly5( out, out_stride, cfg -> phase, m ); break;
        default :
            butterflyN( out, out_stride, cfg -> phase, m, p, cfg -> n );
    }
}

/**
 * @brief Factorizes data size
 * 
 * @param cfg   FFT configuration variable that will hold the results
 * @return      0 on success, -1 on error (unable to factorize)
 */
int fft_init_factors( t_fft_cfg* cfg )
{
    unsigned int i, p;
    unsigned int n = cfg -> n;
    const unsigned int s = floor( sqrt((double) cfg -> n));


    for(i=0, p=4; n>1; i++) {
        if ( i >= MAX_FACTORS ) {
            fprintf(stderr, "(*error*) Unable to factorize fft size.\n");
            return -1;
        }

        while ( n % p ) {
            switch (p) {
                case (4) : p = 2; break;
                case (2) : p = 3; break;
                default  : p += 2;
            }
            if ( p > s ) p = n;
        }

        n = n / p;
        cfg->factors[i].p = p;
        cfg->factors[i].n = n;

        // printf(" %d | %d \n", p, n);
    }

    return 0;
}

/**
 * @brief Initialize FFT configuration
 * 
 * @param cfg           FFT configuration
 * @param n             Number of points
 * @param direction     FFT direction {FFT_FORWARD, FFT_BACKWARD}
 * @return              0 on success, 1 on error
 */
int fft_init_cfg( t_fft_cfg* cfg, unsigned int n, enum fft_direction direction ){

    cfg -> n = n;
    cfg -> direction = direction;


    if ( ! ( cfg -> phase = (float complex *) malloc( n * sizeof(float complex)) ) ) {
        fprintf(stderr, "(*error*) Unable to allocate temporary memory for fft\n");
        return -1;
    }

    double phase_mult =  ( direction == FFT_BACKWARD ) ? (2 * M_PI) / n : - (2 * M_PI) / n;
    for( unsigned int i=0; i<n; i++) {
        cfg -> phase[i] = cos( i * phase_mult ) + I * sin( i * phase_mult );
    }

    return fft_init_factors( cfg );
}

/**
 * @brief Cleanup FFT configuration
 * 
 * @param cfg   FFT configuration
 * @return      0 on success (always returns 0)
 */
int fft_cleanup_cfg( t_fft_cfg* cfg ){
    free( cfg -> phase );
    cfg -> n = 0;

    return 0;
}

/**
 * @brief Perform a Complex-to-Complex FFT
 * 
 * When transforming a complex dataset f of size {N} the output will be a
 * complex dataset F of the same size. The output data is organized as
 * follows:
 * 
 *   * F[ 0 ] : 0 frequency (DC) component
 *   * F[ 0 < j < N/2 ] : $f = j \times \Delta f$
 *   * F[ N/2 ] : $f = \pm f_{Nyquist}$
 *   * F[ N/2 < j < N ] : $f = (j - N) \times \Delta f$
 * 
 * Where $f_{Nyquist} = \frac{1}{2 \Delta t}$ and
 * $\Delta f = \frac{1}{N \Delta t}$ for a signal sampled at $\Delta t$
 * intervals.
 * 
 * @param cfg   FFT configuration
 * @param in    Input data
 * @param out   Output data
 */
void fft_c2c( t_fft_cfg* cfg, const float complex in[], float complex out[] ) {

    transform( in, 1, out, 1, cfg -> factors, cfg );

}

/**
 * @brief Gets the spectral spacing (dk) of points in the transform
 * 
 * @param n     Number of points in input data
 * @param dx    Cell size in real space
 * @return      Cell size in Fourier space
 */
float fft_dk( const unsigned int n, const float dx) {

    float dk = (2 * M_PI) / ( n * dx );

    return dk;
}

/***************************************************************************************
 1D FFT of real data
****************************************************************************************/

/**
 * @brief Initialize real data FFT configuration
 * 
 * @param rcfg          Real data FFT configuration
 * @param nr            Number of points in input (real) data
 * @param direction     Transform direction {FFT_FORWARD, FFT_BACKWARD}
 * @return              0 on sucess, -1 on error
 */
int fftr_init_cfg( t_fftr_cfg* rcfg, unsigned int nr, enum fft_direction direction )
{
    unsigned int i, n;

    if ( nr % 2 ) {
        fprintf(stderr, "(*error*) real ffts are implemented for even sized arrays only.\n");
        return -1;
    }

    n = nr / 2;

    fft_init_cfg( &rcfg -> cfg, n, direction);

    if ( ! ( rcfg -> phase = (float complex *) malloc( (n/2) * sizeof(float complex)) ) ) {
        fprintf(stderr, "(*error*) Unable to allocate temporary memory for fft\n");
        return -1;
    }

    const double p =  ( direction == FFT_BACKWARD ) ? M_PI : - M_PI;
    for( i=0; i<n/2; i++) {
        double phase = p * ( (double) (i+1)/n + 0.5 );
        rcfg -> phase[i] = cos( phase ) + I * sin( phase );
    }

    return 0;
}

/**
 * @brief Cleanup real data FFT configuration
 * 
 * @param rcfg  Real data FFT configuration
 * @return      0 on success (always returns 0)
 */
int fftr_cleanup_cfg( t_fftr_cfg* rcfg ){
    fft_cleanup_cfg( & rcfg -> cfg );
    free( rcfg -> phase );

    return 0;
}

/**
 * @brief Performs a real-to-Complex FFT
 * 
 * A forward Fourier transform of real data will be a conjugate even
 * sequence, meaning that:
 * $\tilde{F}(-k) = \tilde{F}^{*}(k)$
 * This also implies that the imaginary parts of $\tilde{F}(0)$ and
 * $\tilde{F}(\pm f_K)$ must be zero.
 * 
 * Data is stored using a "Complex conjugate storage organization" meaning
 * that for a real dataset with n points we use n/2+1 complex values for
 * the transform, storing only positive frequencies:
 * 
 *   * F[ 0 ] : 0 frequency (DC) component
 *   * F[ 0 < j <= N/2 ] : $f = j \times \Delta f$
 * 
 * If negative frequency component values are required they can be obtained
 * by taking the complex conjugate of the positive counterpart:
 * F[-j] = conj(F[j])
 * 
 * @param rcfg  Real data FFT configuration (must have direction FFT_FORWARD)
 * @param in    Real data input
 * @param out   Complex data output
 */
void fftr_r2c( t_fftr_cfg* rcfg, const float in[], float complex out[] )
{
    if ( rcfg -> cfg.direction != FFT_FORWARD ) {
        fprintf(stderr, "(*error*) configuration for fftr_r2c must have direction = FFT_FORWARD\n");
        return;
    }

    const unsigned int n = rcfg -> cfg.n;
    float complex *buffer = (float complex *) malloc( n * sizeof(float complex) );

    fft_c2c( &rcfg -> cfg, (const float complex *) in, buffer );

    // The first (DC) and last (Nyquist frequency) points are purely real
    // and are calculated from the DC points of the complex transform
    out[0] = creal(buffer[0]) + cimag(buffer[0]);
    out[n] = creal(buffer[0]) - cimag(buffer[0]);

    unsigned int i;
    for( i = 1; i <= n/2; i++) {
        float complex z0, z1, s0, s1;

        z0 = buffer[i];
        z1 = conjf( buffer[n-i] );

        s0 = z0 + z1;
        s1 = (z0 - z1) * rcfg -> phase[i-1];

        out[i]   = 0.5f * ( s0 + s1 );
        out[n-i] = 0.5f * conjf( s0 - s1 );
    }

    free( buffer );

}

/**
 * @brief Perform a Complex to Real FFT
 * 
 * Input data is expected to be organized according to the output of the
 * `fftr_r2c()` routine
 * 
 * @param rcfg  Real data FFT configuration (must have direction FFT_BACKWARD)
 * @param in    Complex data input
 * @param out   Real data output
 */
void fftr_c2r( t_fftr_cfg* rcfg, const float complex in[], float out[] ) {

    if ( rcfg -> cfg.direction != FFT_BACKWARD ) {
        fprintf(stderr, "(*error*) configuration for fftr_c2r must have direction = FFT_BACKWARD\n");
        return;
    }

    const unsigned int n = rcfg -> cfg.n;
    float complex *buffer = (float complex *) malloc( n * sizeof(float complex) );

    buffer[0] = ( creal(in[0]) + creal(in[n]) ) +
            I * ( creal(in[0]) - creal(in[n]) );

    unsigned int i;

    for( i = 1; i <= n/2; i++) {
        float complex z0, z1, s0, s1;
        z0 = in[i];
        z1 = conjf( in[n-i]);

        s0 = z0 + z1;
        s1 = (z0 - z1) * rcfg -> phase[i-1];

        buffer[i  ] =      ( s0 + s1 );
        buffer[n-i] = conjf( s0 - s1 );
    }


    const float norm = 1.0 / (2*n);
    for( i = 0; i < n; i++) buffer[i] *= norm;


    fft_c2c( &rcfg -> cfg, buffer, (float complex *) out);

    free( buffer );
}

/***************************************************************************************
 2D FFT of real data
****************************************************************************************/

/**
 * @brief Initialize 2D real data FFT configuration
 * 
 * @param cfg           2D real data FFT configuration
 * @param nx            Number of x points in input (real) data
 * @param ny            Number of y points in input (real) data
 * @param stridey       Data stride along y
 * @param direction     Transform direction {FFT_FORWARD, FFT_BACKWARD}
 * @return              0 on sucess, -1 on error
 */
int fftr2d_init_cfg( t_fftr2d_cfg* cfg, unsigned int nx, unsigned int ny,
    unsigned int stridey, enum fft_direction direction ) {

    int ierr;

    cfg -> nx = nx;
    cfg -> ny = ny;

    if ( stridey > 0 ) {
      cfg -> stridey = stridey;
    } else {
      cfg -> stridey = nx;
    }

    if ((ierr = fftr_init_cfg( &cfg -> cfgx, nx, direction ))) {
      return(ierr);
    }

    if ((ierr = fft_init_cfg( &cfg -> cfgy, ny, direction ))) {
      return(ierr);
    }

    return(0);
}

/**
 * @brief Cleanup 2D real data FFT configuration
 * 
 * @param cfg   Real data FFT configuration
 * @return      0 on success (always returns 0)
 */
int fftr2d_cleanup_cfg( t_fftr2d_cfg* cfg ){
  fftr_cleanup_cfg( &cfg->cfgx);
  fft_cleanup_cfg( &cfg->cfgy);
  cfg -> nx = cfg -> ny = 0;

  return(0);
}

/**
 * @brief Performs a 2D real-to-Complex FFT
 * 
 * The 2D real to complex transform begins by doing a (1D) real to complex
 * transform of every line along the first dimension, transposing the data
 * and then doing a complex to complex transform of every line of the
 * resulting dataset.
 * 
 * For an input array f(x,y) of dimensions [Nx,Ny] the ouput data will be
 * a transposed array F(ky,kx) of dimensions [Ny, Nx/2+1] organized as:
 * 
 *   * F[ \*, kx ] : $f_x = kx \times \Delta f_x$
 *   * F[ 0, \* ]              : $f_y = 0$
 *   * F[ 0 < ky < Ny/2, \* ]  : $f_y = ky \times \Delta f_y$
 *   * F[ Ny/2, \* ]           : $f_y = \pm f_{Ny}$
 *   * F[ Ny/2 < ky < Ny, \* ] : $f_y = (ky - Ny) \times \Delta f_y$
 *
 * @see fft_c2c()
 * @see fftr_r2c()
 * 
 * @param cfg   Real data FFT configuration (must have direction FFT_FORWARD)
 * @param in    Real data input
 * @param out   Complex data output
 */
void fftr2d_r2c( t_fftr2d_cfg* cfg, const float* in, float complex * out ) {
  unsigned int i,j;

  unsigned const r2c_size = cfg -> nx/2 + 1;

  // Do a real to complex tranform for all the lines and Transpose
  // the result
  for( j = 0; j < cfg -> ny; j++) {
    float complex tmp[ r2c_size ];
    fftr_r2c( &cfg->cfgx, &in[ j* cfg->stridey ], tmp);
    for( i = 0; i < r2c_size; i++ ){
      out[ j + i*cfg -> ny ] = tmp[i];
    }
  }

  // Do an "in-place" complex to complex transform of all the lines
  for( j = 0; j < r2c_size; j++ ) {
    float complex tmp[ cfg -> ny ];
    fft_c2c( &cfg->cfgy, &out[ j* cfg->ny ], tmp);
    for( i = 0; i < cfg -> ny; i++ ) {
      out[ i + j* cfg->ny ] = tmp[i];
    }
  }
}

/**
 * @brief Perform a 2D Complex to Real FFT
 * 
 * Input data is expected to be organized according to the output of the
 * `fftr2d_r2c()` routine
 * 
 * @see fftr2d_r2c()
 * 
 * @param cfg   Real data FFT configuration (must have direction FFT_BACKWARD)
 * @param in    Complex data input
 * @param out   Real data output
 */
void fftr2d_c2r( t_fftr2d_cfg* cfg, const float complex* in, float * out ) {
  unsigned int i,j;

  unsigned const r2c_size = cfg -> nx/2 + 1;

  // Since the size of out[] is less than the size of in[]
  // we cannot use it for the intermediate storage
  float complex tmp[ cfg -> ny * r2c_size ];

  // Normalization for the y back transform; the x back transform is
  // already normalized
  const float norm = 1.0/cfg -> ny;

  // Do a complex to complex transform of all the lines and transpose the
  // result
  for( j = 0; j < r2c_size; j++) {
    float complex tmp_line[ cfg -> ny ];
    fft_c2c( &cfg->cfgy, &in[ j* cfg -> ny ], tmp_line);
    for( i = 0; i < cfg -> ny; i++ ){
      tmp[j + i*r2c_size ] = tmp_line[i] * norm;
    }
  }

  // Do a complex to real transform of all the lines
  for( j = 0; j < cfg -> ny; j++ ) {
    fftr_c2r( &cfg->cfgx, &tmp[ j* r2c_size ], &out[ j * cfg -> stridey ]);
  }
}
