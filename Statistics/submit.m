clc; clear; close all;
addpath('~/Desktop/NI_GIT/Neuroinformatics-M25/Statistics/')
load sampleEEGdata.mat

%% Parameters
chan2use = 'fcz';
min_freq = 3;
max_freq = 30;
num_freq = 20;

chanidx = find(strcmpi({EEG.chanlocs.labels}, chan2use));
freqs   = logspace(log10(min_freq), log10(max_freq), num_freq);

%% Compute time-frequency decomposition
for ev = 1:size(EEG.data,3)
    [cfs(:,:,ev),~] = morletWaveletTransform(EEG.data(chanidx,:,ev), EEG.srate, freqs, 6, 2);
end

power_all = abs(cfs).^2;  % power

%% Trim time
time_s = dsearchn(EEG.times', -500);
time_e = dsearchn(EEG.times', 1200);
eegpower = power_all(:, time_s:time_e, :);
tftimes = EEG.times(time_s:time_e);
nTimepoints = numel(tftimes);

%% Permutation parameters
voxel_pval = 0.01; 
cluster_pval = 0.05;
n_permutes = 2000;
base_idx = [dsearchn(tftimes', -500), dsearchn(tftimes', -100)];

%% Compute baseline power for normalization (will be used in permutations too)
base_power = mean(eegpower(:, base_idx(1):base_idx(2), :), 2); % freqs x 1 x trials

%% Real mean power (baseline normalized in dB, then averaged over trials)
eegpower_dB = 10 * log10(bsxfun(@rdivide, eegpower, base_power)); % dB normalization
realmean = mean(eegpower_dB, 3); % average across trials (dimension 3)

%% Create permutation null distribution
% KEY FIX: Perform circular shift on RAW power data, THEN normalize
permuted_vals = zeros(n_permutes, size(realmean,1), size(realmean,2));

for permi = 1:n_permutes
    % Generate random cutpoint (avoid edges near baseline period)
    cutpoint = randi([2, nTimepoints - diff(base_idx) - 2], 1);
    
    % Circularly shift RAW power data (before normalization)
    shifted_power = eegpower(:, [cutpoint:end, 1:cutpoint-1], :);
    
    % Now compute baseline from the SHIFTED data
    shifted_base_power = mean(shifted_power(:, base_idx(1):base_idx(2), :), 2);
    
    % Normalize the shifted data
    shifted_power_dB = 10 * log10(bsxfun(@rdivide, shifted_power, shifted_base_power));
    
    % Average across trials to get permuted mean
    permuted_vals(permi,:,:) = mean(shifted_power_dB, 3);
end

%% Compute z-map
perm_mean = squeeze(mean(permuted_vals, 1));
perm_std  = squeeze(std(permuted_vals, [], 1));
zmap = (realmean - perm_mean) ./ perm_std;

%% Threshold by voxel p-value (two-tailed test)
voxel_thresh = norminv(1 - voxel_pval/2); % Two-tailed: divide alpha by 2
threshmean = realmean;
threshmean(abs(zmap) < voxel_thresh) = 0;

%% Cluster-based permutation correction
% Find clusters in real data
real_clusters = bwconncomp(abs(zmap) >= voxel_thresh);
real_cluster_masses = zeros(1, real_clusters.NumObjects);

for ci = 1:real_clusters.NumObjects
    real_cluster_masses(ci) = sum(abs(zmap(real_clusters.PixelIdxList{ci})));
end

% Build null distribution of maximum cluster masses
max_cluster_masses = zeros(n_permutes, 1);

fprintf('Computing cluster-based correction...\n');
for permi = 1:n_permutes
    if mod(permi, 500) == 0
        fprintf('  Permutation %d/%d\n', permi, n_permutes);
    end
    
    % Get permuted z-map
    perm_zmap = squeeze(permuted_vals(permi,:,:));
    perm_zmap = (perm_zmap - perm_mean) ./ perm_std;
    
    % Find clusters in permuted data
    perm_clusters = bwconncomp(abs(perm_zmap) >= voxel_thresh);
    
    if perm_clusters.NumObjects > 0
        perm_cluster_masses = zeros(1, perm_clusters.NumObjects);
        for ci = 1:perm_clusters.NumObjects
            perm_cluster_masses(ci) = sum(abs(perm_zmap(perm_clusters.PixelIdxList{ci})));
        end
        max_cluster_masses(permi) = max(perm_cluster_masses);
    else
        max_cluster_masses(permi) = 0;
    end
end

% Determine significant clusters
cluster_threshold = prctile(max_cluster_masses, 100 * (1 - cluster_pval));
significant_clusters = real_cluster_masses > cluster_threshold;

% Create cluster-corrected map
clustermean = zeros(size(realmean));
for ci = 1:real_clusters.NumObjects
    if significant_clusters(ci)
        clustermean(real_clusters.PixelIdxList{ci}) = realmean(real_clusters.PixelIdxList{ci});
    end
end

%% Print Z-score info
fprintf('\n=== Z-score Statistics ===\n');
fprintf('Mean Z: %.3f\n', mean(zmap(:)));
fprintf('Max Z:  %.3f\n', max(zmap(:)));
fprintf('Min Z:  %.3f\n', min(zmap(:)));
fprintf('Voxel threshold (|Z| > %.2f, p=%.3f, two-tailed)\n', voxel_thresh, voxel_pval);

% Find and print significant time-frequency points
[signif_f, signif_t] = find(abs(zmap) >= voxel_thresh);
fprintf('\nSignificant points (|Z| >= %.2f): %d total\n', voxel_thresh, numel(signif_f));

if numel(signif_f) > 0
    fprintf('\nTop significant points:\n');
    fprintf('Frequency (Hz)\tTime (ms)\tZ-value\t\tPower (dB)\n');
    fprintf('-----------------------------------------------------------\n');
    for i = 1:min(numel(signif_f), 20) % print top 20 
        fprintf('%.2f\t\t%.1f\t\t%.3f\t\t%.3f\n', ...
            freqs(signif_f(i)), tftimes(signif_t(i)), ...
            zmap(signif_f(i), signif_t(i)), ...
            realmean(signif_f(i), signif_t(i)));
    end
end

%% Print cluster correction results
fprintf('\n=== Cluster-Based Correction ===\n');
fprintf('Number of clusters found: %d\n', real_clusters.NumObjects);
fprintf('Cluster mass threshold (p=%.3f): %.2f\n', cluster_pval, cluster_threshold);
fprintf('Significant clusters: %d\n', sum(significant_clusters));

if sum(significant_clusters) > 0
    fprintf('\nSignificant cluster details:\n');
    for ci = 1:real_clusters.NumObjects
        if significant_clusters(ci)
            fprintf('  Cluster %d: mass=%.2f, size=%d voxels\n', ...
                ci, real_cluster_masses(ci), length(real_clusters.PixelIdxList{ci}));
        end
    end
end

%% Plotting
% 1. Observed Power (dB)
figure('Position', [100, 100, 1400, 900]);
subplot(2,2,1);
contourf(tftimes, freqs, realmean, 40, 'linecolor', 'none');
title('Observed Power (dB)');
xlabel('Time (ms)'); ylabel('Frequency (Hz)');
colorbar;
set(gca, 'yscale', 'log');
set(gca, 'ytick', [3 5 8 13 20 30]);
colormap(gca, 'jet');
caxis([-3 3]);
hold on;
xline(0, 'k--', 'LineWidth', 1.5); % Stimulus onset

% 2. Z-map
subplot(2,2,2);
contourf(tftimes, freqs, zmap, 40, 'linecolor', 'none');
title('Z-map (Observed vs Permuted)');
xlabel('Time (ms)'); ylabel('Frequency (Hz)');
colorbar;
set(gca, 'yscale', 'log');
set(gca, 'ytick', [3 5 8 13 20 30]);
colormap(gca, 'jet');
caxis([-5 5]);
hold on;
xline(0, 'k--', 'LineWidth', 1.5);

% 3. Voxel-corrected Thresholded Map
subplot(2,2,3);
contourf(tftimes, freqs, threshmean, 40, 'linecolor', 'none');
title(sprintf('Voxel-corrected (|Z| > %.2f, p=%.3f)', voxel_thresh, voxel_pval));
xlabel('Time (ms)'); ylabel('Frequency (Hz)');
colorbar;
set(gca, 'yscale', 'log');
set(gca, 'ytick', [3 5 8 13 20 30]);
colormap(gca, 'jet');
caxis([-3 3]);
hold on;
xline(0, 'k--', 'LineWidth', 1.5);

% 4. Cluster-corrected Map
subplot(2,2,4);
contourf(tftimes, freqs, clustermean, 40, 'linecolor', 'none');
title(sprintf('Cluster-corrected (p=%.3f)', cluster_pval));
xlabel('Time (ms)'); ylabel('Frequency (Hz)');
colorbar;
set(gca, 'yscale', 'log');
set(gca, 'ytick', [3 5 8 13 20 30]);
colormap(gca, 'jet');
caxis([-3 3]);
hold on;
xline(0, 'k--', 'LineWidth', 1.5);

sgtitle(sprintf('Time-Frequency Analysis: Channel %s', upper(chan2use)), 'FontSize', 14, 'FontWeight', 'bold');
% 
% - Computes Morlet wavelet transform (3-30 Hz, 20 frequencies, 6 cycles)
% - Trims time window (-500 to 1200 ms)
% - Performs 2000 circular-shift permutations on raw power
% - Applies baseline normalization (-500 to -100 ms) to each permutation
% - Builds null distribution (mean and std)
% - Computes Z-map comparing observed to null
% - Applies voxel-level threshold (p < 0.01, two-tailed)
% - Implements cluster-based correction (p < 0.05)
% - Identifies significant clusters based on cluster mass statistic
