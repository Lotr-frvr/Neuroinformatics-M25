# EEG Time–Frequency Analysis with Morlet Wavelet Convolution

### 1. Morlet wavelet 

A Morlet wavelet is a complex sinusoid multiplied by a Gaussian window:

$$
	ext{wavelet}(t) = e^{2\pi i f t} \cdot e^{-t^{2}/(2\sigma^{2})}
$$

- The complex sinusoid, $e^{2\pi i f t}$, oscillates at center frequency $f$.
- The Gaussian window, $e^{-t^{2}/(2\sigma^{2})}$, localises the sinusoid in time.
- The resulting wavelet is a time-limited oscillation useful for extracting frequency-specific activity while retaining temporal information.

### 2. Why use complex wavelets?

- The real part of the complex convolution approximates a bandpass-filtered version of the signal at frequency $f$.
- The imaginary part is the 90° phase-shifted component. 
- Magnitude, `abs(convolution)`, gives instantaneous amplitude/power at $f$; `angle(convolution)` gives instantaneous phase as usual.

### 3. Convolution = Matching-Wavelet filtering

Convolving a wavelet with EEG data acts as a matched filter tuned to the wavelet's center frequency and width. At each time point the convolution quantifies similarity between the data and the time-limited oscillatory template.


### 4. FFT-based convolution

Time-domain convolution is equivalent to multiplication in the frequency domain. For long signals, FFT-based convolution is far more efficient:

$$
	ext{convolution} = \text{ifft}\big(\text{fft}(\text{wavelet}) \times \text{fft}(\text{signal})\big)
$$


### `lec_12.m`

- Builds single Morlet wavelets at selected frequencies.
- Visualises the multiplicative product (sine × Gaussian) and shows the real/imaginary components.
- Creates many wavelets to show  changing width and temporal spread.
- Compares wavelet responses with simple bandpass-filtered signals to highlight differences.

### `lec_13.m` 

- Demonstrates Euler's formula and complex-valued wavelets in 3D (real vs imaginary vs time).
- Shows phase extraction and the relationships between filtered signal, power, and phase:
	- Filtered signal: `real(convolution)`
	- Power: `abs(convolution).^2`
	- Phase: `angle(convolution)`
- Introduces frequency scaling, wavelet width metrics (e.g., FWHM), and compares fixed vs adaptive cycles.


##### note: Wavelet cycles determine how many oscillations fit inside the wavelet’s Gaussian window, controlling its time–frequency resolution. Cycle scaling increases the number of cycles with frequency, ensuring low frequencies use longer wavelets and high frequencies use shorter, more precise wavelets