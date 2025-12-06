
%% FFT of EEG Signal 

clear; close all; clc;

%% Load data
load sampleEEGdata
signal = double(EEG.data(10,:,5));
time_eeg = EEG.times;
srate = EEG.srate;

%% FFT
N  = length(signal);
fft_sig = fft(signal)/N;
hz = linspace(0, srate/2, floor(N/2)+1);

%% Plot
figure
plot(hz, abs(fft_sig(1:length(hz))*2), 'linew',2)
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('FFT of EEG signal (single channel)')

%% Gaussian Convolution 



%% Create Gaussian kernel
time = -1:1/srate:1;
s = 5/(2*pi*30);
gaussian = exp((-time.^2)/(2*s^2));
gaussian = gaussian / sum(gaussian);

%% Plot Gaussian
figure
plot(time, gaussian)
title('Gaussian kernel')
xlabel('Time (s)')
ylabel('Amplitude')

%% Convolution
conv_result = conv(signal, gaussian, 'same');

%% Plot raw vs smoothed EEG
figure
subplot(211)
plot(time_eeg, signal)
title('Raw EEG')

subplot(212)
plot(time_eeg, conv_result,'r','linew',2)
title('Gaussian-smoothed EEG')
xlabel('Time (ms)')
