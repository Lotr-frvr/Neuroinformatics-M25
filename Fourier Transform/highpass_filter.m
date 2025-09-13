clc
clear
close all

% ===== Load EEG data =====
load('sampleEEGdata.mat');
signal = EEG.data(10,:,5); % 10th channel, 5th trial
fs = EEG.srate;
t = (0:length(signal)-1) / fs;

%% ===== Butterworth High-Pass Filter =====
% Filter parameters
cutoff = 1;       % cutoff frequency in Hz (typical EEG high-pass = 0.1–1 Hz)
order = 4;        % filter order (higher = steeper slope)

% Normalized cutoff (Nyquist = fs/2)
Wn = cutoff / (fs/2);

% Design Butterworth high-pass
[b,a] = butter(order, Wn, 'high');

% Apply filter (zero-phase using filtfilt avoids delay)
signal_hp = filtfilt(b, a, signal);

%% ===== Compare with Original =====
figure;
subplot(2,1,1);
plot(t, signal, 'b'); hold on;
plot(t, signal_hp, 'r');
xlim([0 5]); % zoom into first 5 seconds
xlabel('Time (s)'); ylabel('Amplitude (µV)');
legend('Raw EEG','High-pass filtered');
title('Butterworth High-Pass Filtering (Time Domain)');

%% ===== Frequency Response of the Filter =====
N = 1024;
[H,f] = freqz(b,a,N,fs);

subplot(2,1,2);
plot(f, abs(H), 'k','LineWidth',1.5);
xlim([0 10]); ylim([0 1.1]);
xlabel('Frequency (Hz)'); ylabel('Gain');
title('Butterworth High-Pass Filter Frequency Response');
grid on;
