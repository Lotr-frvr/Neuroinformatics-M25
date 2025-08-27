### Pre-Processing 

In the reviewed code , the following order is observed:

- Load EEG data
- Apply Notch Filter (for removing 60Hz frequency band)
- Plot PSD
- Re-reference using Common Average Reference
- Epoching
- Baseline Correction
- Artifact Rejection: first trials then channels

There are some small issues with this order of preprocessing

- A lot of processing is done on artifacts that are finally being rejected, though it leads to the intended results, this leads to some redundancy.
- When plotting the PSD, you can identify the bad channels if they do exist. However, the bad channels were removed after other processing stpes.

#### Recommended Pipeline 

A more efficient and widely recommended preprocessing pipeline is:

1. Load EEG data
2. Filtering: (bandpass: e.g., 0.5–40 Hz; notch if needed)
3. Re-reference (average, mastoid, or linked ears)
4. Bad channel removal (variance check, interpolation if needed)
5. Plot PSD (after bad channel removal, can see if there are any remaining anomalies)
6. Epoching (cut into trials based on events)
7. Baseline correction (using pre-stimulus window)
8. Artifact rejection (first trials, then channels if needed)

Key differences:
- Bad channel removal is done before epoching and plotting, reducing redundant processing and avoiding plotting with noisy channels.
- Filtering is broad (bandpass + notch), not just notch. A bandpass keeps only the frequency of interest and removes all the noise.

