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



%% Baseline normalization per trial (in dB)
base_power = mean(eegpower(:, base_idx(1):base_idx(2), :), 2);
eegpower_dB = 10 * log10(bsxfun(@rdivide, eegpower, base_power)); % @rdivide 



%% Real mean power (averaged over trials) 
realmean = mean(eegpower_dB, 3); % 3


%% Create permutation null distribution
permuted_vals = zeros(n_permutes, size(realmean,1), size(realmean,2));

for permi = 1:n_permutes
    cutpoint = randi([2, nTimepoints - diff(base_idx) - 2], 1);
    shifted_data = eegpower_dB(:, [cutpoint:end, 1:cutpoint-1], :); % shifted
    permuted_vals(permi,:,:) = mean(shifted_data, 3); % average across trials
end

%% Compute z-map
perm_mean = squeeze(mean(permuted_vals, 1));
perm_std  = squeeze(std(permuted_vals, [], 1));
zmap = (realmean - perm_mean) ./ perm_std;


%% Threshold by voxel p-value
voxel_thresh = norminv(1 - voxel_pval);
threshmean = realmean;
threshmean(abs(zmap) < voxel_thresh) = 0;



%% Pritn Z score info
fprintf('Mean Z: %.3f\n', mean(zmap(:)));
fprintf('Max Z:  %.3f\n', max(zmap(:)));
fprintf('Min Z:  %.3f\n', min(zmap(:)));
% |Z|>=
fprintf('Voxel threshold (|Z| > %.2f, p=%.3f)\n', voxel_thresh, voxel_pval);

% Find and print significant time-frequency points
[signif_f, signif_t] = find(abs(zmap) >= voxel_thresh);
% |Z|>= 
fprintf('\nSignificant points (|Z| >= %.2f): %d total\n', voxel_thresh, numel(signif_f));
% just printing 
fprintf('Frequency (Hz)\tTime (ms)\tZ-value\n');
for i = 1:min(numel(signif_f), 20) % print top 20 
    fprintf('%.2f\t\t%.1f\t\t%.3f\n', freqs(signif_f(i)), tftimes(signif_t(i)), zmap(signif_f(i), signif_t(i)));
end



%% Plotting
figure;
% Observed Power (dB)
figure;
contourf(tftimes, freqs, realmean, 40, 'linecolor', 'none');
title('Observed Power (dB)');
xlabel('Time (ms)'); ylabel('Frequency (Hz)');
colorbar;
set(gca, 'yscale', 'log');
colormap jet;

%  Z-map
figure;
contourf(tftimes, freqs, zmap, 40, 'linecolor', 'none');
title('Z-map (Observed vs Permuted)');
xlabel('Time (ms)'); ylabel('Frequency (Hz)');
colorbar;
set(gca, 'yscale', 'log');
colormap jet;

% Thresholded Map
figure;
contourf(tftimes, freqs, threshmean, 40, 'linecolor', 'none');
title(['Thresholded Map (|Z| > ' num2str(voxel_thresh,'%.2f') ', p=' num2str(voxel_pval) ')']);
xlabel('Time (ms)'); ylabel('Frequency (Hz)');
colorbar;
set(gca, 'yscale', 'log');
colormap jet;