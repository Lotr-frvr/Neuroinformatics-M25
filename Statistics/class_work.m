clc; close all;clear all;
load sampleEEGdata.mat

% Input parameter 
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

% Compute time-frequency decomposition
for fi = 1:num_freq
    % Create Morlet wavelet
    s = wave_num / (2 * pi * freqs(fi)); % standard deviation
    wavetime = -2:1/EEG.srate:2; % wavelet time vector
    morlet_wave = exp(2*1i*pi*freqs(fi).*wavetime) .* exp(-wavetime.^2./(2*s^2));
    
    % FFT parameters
    n_wavelet = length(morlet_wave);
    n_data = EEG.pnts;
    n_conv = n_wavelet + n_data - 1;
    half_wave = floor(n_wavelet/2);
    
    % FFT of wavelet
    fft_wave = fft(morlet_wave, n_conv);
    
    % Convolve with each trial
    for ti = 1:EEG.trials
        % FFT of data
        fft_data = fft(data(:, ti), n_conv);
        
        % Convolution in frequency domain
        convolution = ifft(fft_wave .* fft_data);
        convolution = convolution(half_wave+1:end-half_wave);
        
        % Store power
        power_all(fi, :, ti) = abs(convolution).^2;
    end
end

fprintf('TFR computed: %d freqs x %d timepoints x %d trials\n', ...
        size(power_all, 1), size(power_all, 2), size(power_all, 3));

%% Initialize null hypothesis and create Z-score

% Define baseline period (e.g., -500 to -200 ms before stimulus)
baseline_time = [-500 -200]; % in ms
baseline_idx = find(EEG.times >= baseline_time(1) & EEG.times <= baseline_time(2));

% Compute baseline mean and std for each frequency across all trials
baseline_power = power_all(:, baseline_idx, :); % freqs x baseline_times x trials
baseline_power_mean = mean(baseline_power(:, :), 2); % freqs x 1
baseline_power_std = std(baseline_power(:, :), 0, 2); % freqs x 1

% Convert power to dB
power_db = 10 * log10(power_all);
baseline_db_mean = 10 * log10(baseline_power_mean);
baseline_db_std = 10 * log10(baseline_power_std);

% Compute Z-scores: (power - baseline_mean) / baseline_std
% Z-score normalization for each frequency
z_power = zeros(size(power_all));
for fi = 1:num_freq
    z_power(fi, :, :) = (power_all(fi, :, :) - baseline_power_mean(fi)) / baseline_power_std(fi);
end

fprintf('Z-score normalization completed\n');

%% Visualization
figure('Position', [100 100 1400 800]);

% Plot 1: Mean power across trials (dB)
subplot(2, 3, 1);
mean_power_db = mean(power_db, 3);
imagesc(EEG.times, freqs, mean_power_db);
axis xy;
colorbar;
xlabel('Time (ms)');
ylabel('Frequency (Hz)');
title('Mean Power (dB)');
colormap('jet');

% Plot 2: Mean Z-scored power
subplot(2, 3, 2);
mean_z_power = mean(z_power, 3);
imagesc(EEG.times, freqs, mean_z_power);
axis xy;
colorbar;
xlabel('Time (ms)');
ylabel('Frequency (Hz)');
title('Mean Z-scored Power');
colormap('jet');

% Plot 3: Baseline-corrected power (dB relative to baseline)
subplot(2, 3, 3);
power_db_baseline_corrected = zeros(size(power_db));
for fi = 1:num_freq
    power_db_baseline_corrected(fi, :, :) = power_db(fi, :, :) - baseline_db_mean(fi);
end
mean_power_db_corrected = mean(power_db_baseline_corrected, 3);
imagesc(EEG.times, freqs, mean_power_db_corrected);
axis xy;
colorbar;
xlabel('Time (ms)');
ylabel('Frequency (Hz)');
title('Baseline-corrected Power (dB)');
colormap('jet');
caxis([-3 3]);

% Plot 4: Power at a specific frequency (e.g., 10 Hz)
subplot(2, 3, 4);
[~, freq_idx] = min(abs(freqs - 10)); % Find closest to 10 Hz
plot(EEG.times, squeeze(mean(power_all(freq_idx, :, :), 3)));
xlabel('Time (ms)');
ylabel('Power');
title(sprintf('Power at %.1f Hz', freqs(freq_idx)));
grid on;

% Plot 5: Z-score at same frequency
subplot(2, 3, 5);
plot(EEG.times, squeeze(mean(z_power(freq_idx, :, :), 3)));
xlabel('Time (ms)');
ylabel('Z-score');
title(sprintf('Z-scored Power at %.1f Hz', freqs(freq_idx)));
grid on;
hold on;
yline(0, 'k--', 'LineWidth', 1);
yline(1.96, 'r--', 'p < 0.05');
yline(-1.96, 'r--', 'p < 0.05');

% Plot 6: Frequency spectrum at a specific time (e.g., 300 ms)
subplot(2, 3, 6);
[~, time_idx] = min(abs(EEG.times - 300)); % Find closest to 300 ms
plot(freqs, squeeze(mean(z_power(:, time_idx, :), 3)));
xlabel('Frequency (Hz)');
ylabel('Z-score');
title(sprintf('Z-scored Power at %d ms', round(EEG.times(time_idx))));
grid on;
set(gca, 'XScale', 'log');
hold on;
yline(0, 'k--', 'LineWidth', 1);
yline(1.96, 'r--', 'p < 0.05');

sgtitle(sprintf('Time-Frequency Analysis - Channel: %s', chan2use)); 