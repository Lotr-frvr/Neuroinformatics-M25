# Fourier Analysis, Noise, and Filtering of Time Series and EEG Data



### 1 Sine waves 

A sine wave :

$$
x(t) = A \sin(2\pi f t + \phi)
$$

### 2 DFT and FFT

For a discrete time series $x[n]$ of length $N$, the DFT is:

$$
X[k] = \sum_{n=0}^{N-1} x[n] e^{-i 2\pi k n / N}
$$

The frequency corresponding to bin $k$ is $f_k = k\cdot \frac{srate}{N}$. The power is the squared magnitude of the complex spectrum:

$$
P[k] = |X[k]|^2
$$


### 3 Effect of noise: time vs. frequency domain

Additive noise (especially white noise) distributes energy across frequencies.

In the time domain noise can obscure  structure, but in the frequency domain true oscillatory components appear as localised peaks (narrowband increases in power) on top of a broadband noise floor. Also, spectral methods are often more sensitive for detecting oscillations in noisy data.

- Small noise: time-domain waveform still shows oscillations; spectrum shows sharp peaks.
- Large noise: time-domain oscillations are obscured; spectrum still reveals peaks at true frequencies.

### 4 Convolution and the frequency domain

Discrete convolution is defined as:

$$
y[n] = \sum_{k} x[k] \cdot h[n-k]
$$

convolution in time corresponds to multiplication in frequency. Therefore, designing a kernel $h$ with a desired frequency response is equivalent to designing a filter.

### 5 Gaussian convolution as a low-pass filter

A Gaussian kernel in time is:

$$
g(t) = e^{-t^2 / (2\sigma^2)}
$$

Its Fourier transform is also Gaussian. Convolution with a Gaussian thus performs smoothing (low-pass filtering), attenuating high-frequency noise while preserving slow components.

The Fourier Transform of a Gaussian is another Gaussian, meaning it has energy at all frequencies but with the highest energy (and largest weights) around zero frequency (DC) and decreasing energy as frequency increases.

### 6 High-pass filters: difference kernel and FIR

Difference kernel (high-pass):

$$
h[n] = [-1, \; 1]
$$

This computes a discrete derivative: $y[n] = x[n] - x[n-1]$, removing DC and very slow components. When yo convolve $ h[n] $ with $x[n]$ you get $y[n]$

Finite Impulse Response (FIR) high-pass filters are designed to have controlled magnitude and linear-phase properties. Steps:

- Choose cutoff frequency $f_c$ (Hz).
- Choose filter order $N$ (tradeoff: transition bandwidth vs. ripple).
- Window and truncate an ideal high-pass impulse response, or use library functions (e.g., MATLAB `fir1`).

Applying an FIR filter with zero-phase filtering (e.g., MATLAB `filtfilt`) removes phase distortions and preserves the temporal shape of oscillatory activity while removing slow drifts.

### 7 Why FFT is preferred over a direct DFT

Computing the DFT directly (matrix multiplication) scales as $O(N^2)$, which becomes prohibitively slow for long signals or when repeated on many trials/channels. The FFT reduces this to $O(N \log N)$, permitting fast spectral estimates, real-time applications, and practical permutation/bootstrap methods that require many transforms.

### 8 Time vs frequency domain interpretation

- Time domain: shows waveform morphology and precise timing of events; sensitive to transient, non-stationary features.
- Frequency domain: summarises the energy distribution across frequencies; more sensitive to sustained oscillatory structure and useful for detecting periodicities masked by noise.
