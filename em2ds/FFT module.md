# FFT module

## Frequency values

The discrete Fourier transform of a signal sampled at $\Delta t$ intervals will hold frequencies in the $[-f_{Nyquist},+f_{Nyquist}]$ interval, where $f_{Nyquist}$ is the Nyquist frequency, defined as half the sampling rate or:

$f_{Nyquist} = \frac{1}{2 \Delta t}$

If the signal sample has $N$ points then the discrete transform points will be spaced at $\Delta f$ intervals defined as:

$\Delta f = \frac{1}{N \Delta t}$

or in terms of angular frequency $\Delta \omega$:

$\Delta \omega = \frac{2 \pi}{N \Delta t}$

## Complex to Complex FFT routine

When transforming a complex dataset f of size {N} the output will be a complex dataset F of the same size. The output data is organized as follows:

* F[ 0 ] : 0 frequency (DC) component
* F[ 0 < j < N/2 ] : $f = j \times \Delta f$
* F[ N/2 ] : $f = \pm f_{Nyquist}$
* F[ N/2 < j < N ] : $f = (j - N) \times \Delta f$

The frequency value can be easily calculated using the following code snippet:

```C
// j must be in the range [0, N[
float df = 1.0f / (N * dt);
float f = ((j <= N/2)? j : (j - (int) N))*df;
```

Note that the this returns the positive Nyquist frequency for point N/2.

## Real to Complex FFT routine

When transforming real data, the FFT routines make use of the fact that the result of the forward transform of real data is a conjugate-even sequence. Consequently, we don't need to store the complete output sequence, and for a real input dataset f of size {N}, we require a complex output dataset F of size {N/2 + 1}.

F[0] holds the 0 frequency component of the transform, and F[N/2] holds the component corresponding to the Nyquist frequency. Due to the symmetry property mentioned above, both these values will be purely real. The frequency value of element F[ j ] can then be calculated using:

```C
// j must be in the range [0, N/2]
float df = 1.0f / (N * dt);
float f = j * df;
```

If negative frequency component values are required they can be obtained by taking the complex conjugate of the positive counterpart:

F[-j] = conj(F[j])

## 2D FFT routines

### 2D Real to Complex FFT routines

The 2D real to complex transform begins by doing a (1D) real to complex transform of every line along the first dimension, transposing the data and then doing a complex to complex transform of every line of the resulting dataset. This means the output of a real to complex 2D transform of a real {Nx,Ny} dataset f using fftr2d_r2c is a transposed complex array F of dimensions {Ny, Nx/2+1}. The data is organized as follows:

* F[ \*, kx ] : $f_x = kx \times \Delta f_x$
* F[ 0, \* ]              : $f_y = 0$
* F[ 0 < ky < Ny/2, \* ]  : $f_y = ky \times \Delta f_y$
* F[ Ny/2, \* ]           : $f_y = \pm f_{Ny}$
* F[ Ny/2 < ky < Ny, \* ] : $f_y = (ky - Ny) \times \Delta f_y$

The frequency values of element F[ ky, kx ] can then be easily calculated using:

```C
// kx must be in the range [0, Nx/2]
// ky must be in the range [0, Ny[
float dfx = 1.0f / (Nx * dx);
float dfy = 1.0f / (Ny * dy);

float fx = kx * dfx;
float fy = ((ky <= Ny/2)? ky : ky - Ny)*dfy;
```
