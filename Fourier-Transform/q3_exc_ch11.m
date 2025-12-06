%% Exercise 3 — Generate & Sum Sine Waves, Add Noise, Compute FFT
% Generate a time series by creating and summing sine waves, as in
% figure 11.2B. Use between two and four sine waves, so that the
% individual sine waves are still somewhat visible in the sum. Perform a
% Fourier analysis (you can use the fft function) on the resulting time
% series and plot the power structure. Confirm that your code is correct by
% comparing the frequencies with nonzero power to the frequencies of the
% sine waves that you generated. Now try adding random noise to the
% signal before computing the Fourier transform. First, add a small
% amount of noise so that the sine waves are still visually recognizable.
% Next, add a large amount of noise so that the sine waves are no longer
% visually recognizable in the time domain data. Perform a Fourier
% analysis on the two noisy signals and plot the results. What is the effect
% of a small and a large amount of noise in the power spectrum? Are the
% sine waves with noise easier to detect in the time domain or in the
% frequency domain, or is it equally easy/difficult to detect a sine wave in
% the presence of noise?

clear; close all; clc;

%% Time and sampling rate
srate = 1000;
time  = 0:1/srate:2;

%% Create sine waves 
frex   = [ 5 12 20 ];
amp = [ 3 4 2 ];
phases = 2*pi*rand(1,3);

sine_waves = zeros(length(frex), length(time));

for fi = 1:length(frex)
    sine_waves(fi,:) = amp(fi) * sin(2*pi*frex(fi).*time + phases(fi));
end

%% Sum the waves
signal_clean = sum(sine_waves);

figure
plot(time, signal_clean, 'k')
title('Sum of three sine waves (clean)')
xlabel('Time (s)')
ylabel('Amplitude')

%% FFT of clean signal
N = length(signal_clean);
hz = linspace(0, srate/2, floor(N/2)+1);

fft_clean = fft(signal_clean)/N;

figure
plot(hz, abs(fft_clean(1:length(hz))*2),'linew',2)
xlabel('Frequency (Hz)')
ylabel('Power')
title('Power spectrum (clean)')

%% SMALL noise
noise_small = randn(size(signal_clean))*1;
signal_small = signal_clean + noise_small;

figure
plot(time, signal_small)
title('Signal with SMALL noise')

fft_small = fft(signal_small)/N;

figure
plot(hz, abs(fft_small(1:length(hz))*2),'linew',2)
xlabel('Frequency (Hz)')
ylabel('Power')
title('Power spectrum (small noise)')

%% LARGE noise
noise_large = randn(size(signal_clean))*8;
signal_large = signal_clean + noise_large;

figure
plot(time, signal_large)
title('Signal with LARGE noise')

fft_large = fft(signal_large)/N;

figure
plot(hz, abs(fft_large(1:length(hz))*2),'linew',2)
xlabel('Frequency (Hz)')
ylabel('Power')
title('Power spectrum (large noise)')
