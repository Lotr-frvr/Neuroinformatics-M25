# EEG Time-Frequency Analysis with Permutation Testing


### Baseline Normalization

**Purpose**: Remove frequency-dependent and individual differences in oscillatory power.

**Decibel (dB) Conversion**:
```
Power_dB(f,t) = 10 × log₁₀(Power(f,t) / Power_baseline(f))
```

**Properties of dB**:
- Converts multiplicative changes to additive
- Compresses dynamic range
- Symmetric: doubling = +3 dB, halving = -3 dB
- 0 dB = no change from baseline

**Baseline Period**: Pre-stimulus window (typically -500 to -200 ms) representing resting neural activity.

### Multiple Comparisons Problem

**Issue**: Testing thousands of time-frequency points increases false positive rate.

**Example**:
- Testing 1000 points at α = 0.05
- Expected false positives = 1000 × 0.05 = 50
- Family-wise error rate approaches 100%

**Solution**: Permutation testing with correction.

### Permutation Testing

**Principle**: Under null hypothesis (no effect), temporal relationships are arbitrary.

**Circular Shift Method**:
```
For each permutation:
  1. Choose random cutpoint in time
  2. Circularly shift raw power: [cutpoint:end, 1:cutpoint-1]
  3. Apply baseline normalization to shifted data
  4. Compute mean power
```

**Why circular shift?**
- Preserves autocorrelation structure
- Destroys stimulus-locked timing
- Creates valid null distribution

**Critical**: Must shift RAW power before normalization, not after.

**Null Distribution Construction**:
```
After n_permutes iterations:
  μ_null(f,t) = mean(permuted_values)
  σ_null(f,t) = std(permuted_values)
```

### Z-Score Normalization

**Formula**:
```
Z(f,t) = (Observed(f,t) - μ_null(f,t)) / σ_null(f,t)
```

**Interpretation**:
- Z = 0: Observed equals null expectation
- Z = 2.33: 2.33 standard deviations above null (p < 0.01, one-tailed)
- Z = 2.58: 2.58 standard deviations above null (p < 0.01, two-tailed)
- Negative Z: Power decrease (desynchronization)


### Cluster-Based Correction

**Principle**: Use spatial-temporal clustering to control family-wise error rate.

**Algorithm**:

1. **Threshold Z-map** at voxel level (e.g., |Z| > 2.58)

2. **Identify clusters** of contiguous significant voxels

3. **Compute cluster mass**:
   ```
   Mass = Σ(|Z-values| in cluster)
   ```

4. **Build null distribution**:
   - For each permutation: compute Z-map, threshold, find max cluster mass
   - Create distribution of maximum cluster masses

5. **Statistical test**:
   ```
   Cluster threshold = percentile(max_masses, 100(1-α))
   p_cluster = (# permutations with mass ≥ observed) / n_permutes
   ```

**Advantages**:
- Controls family-wise error rate
- More sensitive than Bonferroni correction
- Accounts for spatial-temporal structure

### Frequency Bands

**Logarithmic Spacing**:
```matlab
freqs = logspace(log10(min_freq), log10(max_freq), num_freq)
```

- Equal spacing in log space

### Event-Related Synchronization/Desynchronization

**ERS (Event-Related Synchronization)**:
- Positive Z-values / positive dB
- Increased oscillatory power
- Neural populations firing synchronously
- Indicates active processing

**ERD (Event-Related Desynchronization)**:
- Negative Z-values / negative dB
- Decreased oscillatory power
- Reduced synchrony, increased neural diversity
- Indicates cortical activation





