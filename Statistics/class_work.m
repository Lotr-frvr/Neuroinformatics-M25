clc; close all; clear all;
load sampleEEGdata.mat

% Input parameters
chan2use = 'fcz';
min_freq = 3;
max_freq = 30;
num_freq = 20;

%% Get TFR using morlet wavelet of wavenumber 6 
% to get the freqs*time*trials matrix of power. Name it power_all

chanidx = find(strcmpi({EEG.chanlocs.labels}, chan2use));
freqs   = logspace(log10(min_freq), log10(max_freq), num_freq);

% Wavelet parameters
wave_num = 6; % number of wavelet cycles
time = EEG.times / 1000; % convert to seconds

% Get data for the channel
data = squeeze(EEG.data(chanidx, :, :)); % time x trials

% Initialize power matrix: freqs x time x trials
power_all = zeros(num_freq, EEG.pnts, EEG.trials);

fprintf('Computing time-frequency decomposition...\n');

% Compute time-frequency decomposition
for fi = 1:num_freq
    if mod(fi, 5) == 0
        fprintf('  Processing frequency %d/%d (%.2f Hz)\n', fi, num_freq, freqs(fi));
    end
    
    % Create Morlet wavelet
    s = wave_num / (2 * pi * freqs(fi)); % standard deviation
    wavetime = -2:1/EEG.srate:2; % wavelet time vector
    morlet_wave = exp(2*1i*pi*freqs(fi).*wavetime) .* exp(-wavetime.^2./(2*s^2));
    
    % FFT parameters
    n_wavelet = length(morlet_wave);
    n_data = EEG.pnts;
    n_conv = n_wavelet + n_data - 1;
    half_wave = floor(n_wavelet/2);
    
    % FFT of wavelet (only compute once per frequency)
    fft_wave = fft(morlet_wave, n_conv);
    
    % Vectorized convolution across all trials
    % FFT of all trials at once
    fft_data = fft(data, n_conv, 1); % FFT along dimension 1 (time)
    
    % Convolution in frequency domain (broadcasting)
    convolution = ifft(fft_wave.' .* fft_data, [], 1);
    
    % Trim edges and extract valid convolution
    convolution = convolution(half_wave+1:end-half_wave, :);
    
    % Store power for all trials
    power_all(fi, :, :) = abs(convolution).^2;
end

fprintf('TFR computed: %d freqs x %d timepoints x %d trials\n', ...
        size(power_all, 1), size(power_all, 2), size(power_all, 3));

%% Initialize null hypothesis and create Z-score

% Define baseline period (e.g., -500 to -200 ms before stimulus)
baseline_time = [-500 -200]; % in ms
baseline_idx = find(EEG.times >= baseline_time(1) & EEG.times <= baseline_time(2));

fprintf('\nComputing baseline statistics...\n');
fprintf('Baseline period: %.0f to %.0f ms (%d timepoints)\n', ...
    baseline_time(1), baseline_time(2), length(baseline_idx));

%% Method 1: Proper dB-based baseline correction
% Extract baseline power for each frequency and trial
baseline_power = power_all(:, baseline_idx, :); % freqs x baseline_times x trials

% Compute mean baseline power for each frequency across time and trials
% Shape: freqs x 1 x 1 (for broadcasting)
baseline_power_mean = mean(baseline_power, [2, 3]); % mean across time and trials

% Convert to dB relative to baseline
% This is the standard approach in time-frequency analysis
power_db_corrected = 10 * log10(bsxfun(@rdivide, power_all, baseline_power_mean));

fprintf('Baseline-corrected power (dB) computed\n');

%% Method 2: Z-score normalization (trial-by-trial)
% For proper z-scoring, we need baseline statistics in dB space
baseline_db = 10 * log10(baseline_power);

% Compute baseline mean and std for each frequency (across time and trials)
baseline_db_mean = mean(baseline_db, [2, 3]); % freqs x 1 x 1
baseline_db_std = std(baseline_db, 0, [2, 3]); % freqs x 1 x 1

% Convert all power to dB first
power_db_all = 10 * log10(power_all);

% Compute Z-scores: (power_dB - baseline_mean_dB) / baseline_std_dB
z_power = bsxfun(@rdivide, ...
    bsxfun(@minus, power_db_all, baseline_db_mean), ...
    baseline_db_std);

fprintf('Z-score normalization completed\n');

%% Alternative Method 3: Percent change from baseline
% Sometimes used in practice
percent_change = 100 * (power_all - baseline_power_mean) ./ baseline_power_mean;

%% Compute mean values for visualization
mean_power_raw = mean(power_all, 3);
mean_power_db_corrected = mean(power_db_corrected, 3);
mean_z_power = mean(z_power, 3);
mean_percent_change = mean(percent_change, 3);

%% Visualization
figure('Position', [50 50 1600 1000]);

% Plot 1: Raw mean power
subplot(3, 3, 1);
imagesc(EEG.times, freqs, mean_power_raw);
axis xy;
colorbar;
xlabel('Time (ms)');
ylabel('Frequency (Hz)');
title('Mean Raw Power');
set(gca, 'YScale', 'log');
set(gca, 'YTick', [3 5 8 13 20 30]);
colormap(gca, 'parula');
hold on;
xline(0, 'w--', 'LineWidth', 1.5);
xline(baseline_time(1), 'r--', 'LineWidth', 1);
xline(baseline_time(2), 'r--', 'LineWidth', 1);

% Plot 2: Baseline-corrected power (dB) - RECOMMENDED
subplot(3, 3, 2);
imagesc(EEG.times, freqs, mean_power_db_corrected);
axis xy;
colorbar;
xlabel('Time (ms)');
ylabel('Frequency (Hz)');
title('Baseline-corrected Power (dB) [STANDARD]');
set(gca, 'YScale', 'log');
set(gca, 'YTick', [3 5 8 13 20 30]);
colormap(gca, 'jet');
caxis([-3 3]);
hold on;
xline(0, 'w--', 'LineWidth', 1.5);
xline(baseline_time(1), 'k--', 'LineWidth', 0.5);
xline(baseline_time(2), 'k--', 'LineWidth', 0.5);

% Plot 3: Z-scored power
subplot(3, 3, 3);
imagesc(EEG.times, freqs, mean_z_power);
axis xy;
colorbar;
xlabel('Time (ms)');
ylabel('Frequency (Hz)');
title('Z-scored Power');
set(gca, 'YScale', 'log');
set(gca, 'YTick', [3 5 8 13 20 30]);
colormap(gca, 'jet');
caxis([-3 3]);
hold on;
xline(0, 'w--', 'LineWidth', 1.5);
xline(baseline_time(1), 'k--', 'LineWidth', 0.5);
xline(baseline_time(2), 'k--', 'LineWidth', 0.5);

% Plot 4: Percent change
subplot(3, 3, 4);
imagesc(EEG.times, freqs, mean_percent_change);
axis xy;
colorbar;
xlabel('Time (ms)');
ylabel('Frequency (Hz)');
title('Percent Change from Baseline');
set(gca, 'YScale', 'log');
set(gca, 'YTick', [3 5 8 13 20 30]);
colormap(gca, 'jet');
caxis([-100 100]);
hold on;
xline(0, 'w--', 'LineWidth', 1.5);

% Plot 5: Power at a specific frequency (e.g., 10 Hz - Alpha)
subplot(3, 3, 5);
[~, freq_idx] = min(abs(freqs - 10)); % Find closest to 10 Hz
plot(EEG.times, squeeze(mean(power_all(freq_idx, :, :), 3)), 'b-', 'LineWidth', 1.5);
xlabel('Time (ms)');
ylabel('Power (raw)');
title(sprintf('Power at %.1f Hz (Alpha)', freqs(freq_idx)));
grid on;
hold on;
xline(0, 'k--', 'LineWidth', 1.5);
patch([baseline_time(1) baseline_time(2) baseline_time(2) baseline_time(1)], ...
    [ylim fliplr(ylim)], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');

% Plot 6: dB-corrected power at same frequency
subplot(3, 3, 6);
plot(EEG.times, squeeze(mean(power_db_corrected(freq_idx, :, :), 3)), 'b-', 'LineWidth', 1.5);
xlabel('Time (ms)');
ylabel('Power (dB)');
title(sprintf('dB Power at %.1f Hz', freqs(freq_idx)));
grid on;
hold on;
xline(0, 'k--', 'LineWidth', 1.5);
yline(0, 'r--', 'LineWidth', 1);
patch([baseline_time(1) baseline_time(2) baseline_time(2) baseline_time(1)], ...
    [ylim fliplr(ylim)], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');

% Plot 7: Z-score at same frequency with significance lines
subplot(3, 3, 7);
plot(EEG.times, squeeze(mean(z_power(freq_idx, :, :), 3)), 'b-', 'LineWidth', 1.5);
xlabel('Time (ms)');
ylabel('Z-score');
title(sprintf('Z-scored Power at %.1f Hz', freqs(freq_idx)));
grid on;
hold on;
xline(0, 'k--', 'LineWidth', 1.5);
yline(0, 'k-', 'LineWidth', 1);
yline(1.96, 'r--', 'p < 0.05', 'LineWidth', 1);
yline(-1.96, 'r--', 'p < 0.05', 'LineWidth', 1);
yline(2.58, 'r:', 'p < 0.01', 'LineWidth', 1);
yline(-2.58, 'r:', 'p < 0.01', 'LineWidth', 1);
patch([baseline_time(1) baseline_time(2) baseline_time(2) baseline_time(1)], ...
    [ylim fliplr(ylim)], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');

% Plot 8: Frequency spectrum at a specific time (e.g., 300 ms post-stimulus)
subplot(3, 3, 8);
[~, time_idx] = min(abs(EEG.times - 300)); % Find closest to 300 ms
plot(freqs, squeeze(mean(z_power(:, time_idx, :), 3)), 'b-o', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('Frequency (Hz)');
ylabel('Z-score');
title(sprintf('Z-scored Power at %d ms', round(EEG.times(time_idx))));
grid on;
set(gca, 'XScale', 'log');
set(gca, 'XTick', [3 5 8 13 20 30]);
hold on;
yline(0, 'k-', 'LineWidth', 1);
yline(1.96, 'r--', 'p < 0.05', 'LineWidth', 1);
yline(-1.96, 'r--', 'LineWidth', 1);

% Plot 9: Single trial example (raw power)
subplot(3, 3, 9);
trial_to_plot = 1;
imagesc(EEG.times, freqs, power_all(:, :, trial_to_plot));
axis xy;
colorbar;
xlabel('Time (ms)');
ylabel('Frequency (Hz)');
title(sprintf('Single Trial (#%d) Raw Power', trial_to_plot));
set(gca, 'YScale', 'log');
set(gca, 'YTick', [3 5 8 13 20 30]);
colormap(gca, 'parula');
hold on;
xline(0, 'w--', 'LineWidth', 1.5);

sgtitle(sprintf('Time-Frequency Analysis - Channel: %s', upper(chan2use)), ...
    'FontSize', 16, 'FontWeight', 'bold');

%% Statistical summary
fprintf('\n=== Statistical Summary ===\n');
fprintf('dB-corrected power:\n');
fprintf('  Mean: %.3f dB, Std: %.3f dB\n', mean(mean_power_db_corrected(:)), std(mean_power_db_corrected(:)));
fprintf('  Range: [%.3f, %.3f] dB\n', min(mean_power_db_corrected(:)), max(mean_power_db_corrected(:)));

fprintf('\nZ-scored power:\n');
fprintf('  Mean: %.3f, Std: %.3f\n', mean(mean_z_power(:)), std(mean_z_power(:)));
fprintf('  Range: [%.3f, %.3f]\n', min(mean_z_power(:)), max(mean_z_power(:)));

% Find peak power in dB
[max_db, max_idx] = max(mean_power_db_corrected(:));
[max_f, max_t] = ind2sub(size(mean_power_db_corrected), max_idx);
fprintf('\nPeak power (dB): %.3f dB at %.2f Hz, %.1f ms\n', ...
    max_db, freqs(max_f), EEG.times(max_t));

% Find peak z-score
[max_z, max_z_idx] = max(mean_z_power(:));
[max_z_f, max_z_t] = ind2sub(size(mean_z_power), max_z_idx);
fprintf('Peak Z-score: %.3f at %.2f Hz, %.1f ms\n', ...
    max_z, freqs(max_z_f), EEG.times(max_z_t));

fprintf('\n=== Analysis Complete ===\n');
% 
% - Computes Morlet wavelet transform (3-30 Hz, 20 frequencies, 6 cycles)
% - Uses vectorized FFT convolution for efficiency
% - Computes baseline statistics (-500 to -200 ms) across time and trials
% - Applies three normalization methods:
%   - dB-corrected power (standard method)
%   - Z-scored power (for statistics)
%   - Percent change from baseline
% - Generates visualization plots:
%   - Raw power
%   - dB-corrected power                                                                                                                                                        
%   - Z-scored power
%   - Percent change
%   - Alpha band (10 Hz) time courses (raw, dB, Z-score)
%   - Frequency spectrum at 300 ms
%   - Single trial example
% - Shows baseline period and significance thresholds
