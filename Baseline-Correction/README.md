# Time–Frequency Analysis & Baseline Normalization

##  Time–Frequency (TF) Representation
Time–frequency analysis reveals which oscillations are present, when they increase or decrease, and how long they last. 
### How TF data is computed

For each channel the EEG is convolved with a complex Morlet wavelet of the form

$$
W(t,f) = e^{2\pi i f t} \cdot e^{-t^2/(2\sigma_f^2)}
$$

This produces a complex-valued time–frequency signal

$$
Z(t,f) = \text{EEG}(t) * W(t,f)
$$

The power at each time–frequency point is

$$
P(t,f) = |Z(t,f)|^2
$$

### TF matrix dimensions

After computing power for every channel, frequency and time point we obtain a 3‑D TF matrix:

```
channels × frequencies × time
```
Each entry represents oscillatory power at a specific time, frequency and electrode.

## Why Baseline Normalization Is Needed

Raw TF power is difficult to compare across frequencies, channels, and participants because:

- Low frequencies naturally exhibit higher power (the 1/f characteristic).
- Channels and participants have different amplitude scales.
- Noise and non-neural artifacts affect absolute power values.

To compare responses meaningfully we normalize power relative to a baseline period where no stimulus-related activity is assumed (for example, −500 to −200 ms).

For each frequency:

```
baseline_power(f) = mean power during baseline window
```

## dB Normalization 

$$
\text{dB}(t,f) = 10 \log_{10}\left(\frac{P(t,f)}{P_{\text{baseline}}(f)}\right)
$$

Interpretation examples:

| dB value | Meaning |
|---:|---|
| 0 dB | no change from baseline |
| +3 dB | ≈ 40% increase |
| −3 dB | ≈ 30% decrease |
| +10 dB | 10× increase |

Advantages of dB normalization:

- Reduces the 1/f bias across frequencies
- Makes maps comparable across channels and subjects
- Stabilizes variance for subsequent statistics

Percent change instead of dB, compute

$$
\%\Delta(t,f) = 100 \cdot \frac{P(t,f) - P_{\text{baseline}}(f)}{P_{\text{baseline}}(f)}
$$

## 4. Topographical Mapping

To visualise spatial patterns, generate scalp maps at selected time points using TF power across channels. In this exercise we typically show:

- Row 1: raw power maps
- Row 2: dB-normalized power maps


## Summary

- The TF matrix captures oscillatory power across channels, frequencies and time.
- Baseline normalization converts raw power into a relative change metric that is interpretable across frequencies and electrodes.
- dB conversion is recommended for stability and interpretability.
- Topographical maps show how oscillatory responses evolve across the scalp and how baseline correction affects interpretation.
