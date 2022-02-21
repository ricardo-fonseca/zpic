/*
 *  fft.h
 *
 */

#ifndef __FFT__
#define __FFT__

#include <complex.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288   ///< pi
#endif

/**
 * @brief Maximum number of factors for decomposing FFT size
 * 
 */
#define MAX_FACTORS 32

/**
 * @brief FFT direction
 * 
 */
enum fft_direction {
    FFT_FORWARD,    ///< Forward transformation
    FFT_BACKWARD    ///< Backward transformation
};

/**
 * @brief Decomposition factors of the FFT size
 * 
 * For each factor, the remaining size of the FFT will be p*n; the last
 * factor will always have n = 1.
 */
typedef struct {
	unsigned int p; ///< p factor (FFT size will be p * n)
    unsigned int n; ///< n factor (FFT size will be p * n)
} t_fft_factors;

/**
 * @brief FFT Configuration
 * 
 */
typedef struct {
	unsigned int n;                     ///< Grid size
    enum fft_direction direction;       ///< Tranform direction {FFT_FORWARD, FFT_BACKWARD}
    t_fft_factors factors[MAX_FACTORS]; ///< Factorization of Grid size
    float complex *phase;               ///< Phase factors for FFT calculation
} t_fft_cfg;

/**
 * @brief Real data FFT Configuration
 * 
 */
typedef struct {
	t_fft_cfg cfg;          ///< Underlying Complex to Complex FFT configuration
    float complex *phase;   ///< Phase factor for R2C / C2R transforms
} t_fftr_cfg;

/**********************************************************************************
* Complex data FFT routines
***********************************************************************************/
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
void fft_c2c( t_fft_cfg* cfg, const float complex in[], float complex out[] );

/**
 * @brief Initialize FFT configuration
 * 
 * @param cfg           FFT configuration
 * @param n             Number of points
 * @param direction     FFT direction {FFT_FORWARD, FFT_BACKWARD}
 * @return              0 on success, 1 on error
 */
int fft_init_cfg( t_fft_cfg* cfg, unsigned int n, enum fft_direction direction );

/**
 * @brief Cleanup FFT configuration
 * 
 * @param cfg   FFT configuration
 * @return      0 on success (always returns 0)
 */
int fft_cleanup_cfg( t_fft_cfg* cfg );

/**
 * @brief Gets the spectral spacing (dk) of points in the transform
 * 
 * @param n     Number of points in input data
 * @param dx    Cell size in real space
 * @return      Cell size in Fourier space
 */
float fft_dk( const unsigned int n, const float dx);

/**********************************************************************************
* Real data FFT routines
***********************************************************************************/

/**
 * @brief Initialize real data FFT configuration
 * 
 * @param rcfg          Real data FFT configuration
 * @param nr            Number of points in input (real) data
 * @param direction     Transform direction {FFT_FORWARD, FFT_BACKWARD}
 * @return              0 on sucess, -1 on error
 */
int fftr_init_cfg( t_fftr_cfg* rcfg, unsigned int nr, enum fft_direction direction );

/**
 * @brief Cleanup real data FFT configuration
 * 
 * @param rcfg  Real data FFT configuration
 * @return      0 on success (always returns 0)
 */
int fftr_cleanup_cfg( t_fftr_cfg* rcfg );

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
 * that for a real dataset with N points we use N/2+1 complex values for
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
void fftr_r2c( t_fftr_cfg* rcfg, const float in[], float complex out[] );

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
void fftr_c2r( t_fftr_cfg* rcfg, const float complex in[], float out[] );

#endif
