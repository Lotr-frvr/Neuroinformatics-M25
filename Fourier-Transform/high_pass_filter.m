%% high_pass_filter
% -  A: simple difference kernel (for demonstration)
% -  B: FIR high-pass designed with fir1 and applied with filtfilt (recommended)
%

clear; close all; clc;

%% Load example EEG (same as in your other files)
load sampleEEGdata;        % provides EEG struct
signal = double(EEG.data(10,:,5));  % 10th channel, trial 5 (as before)
srate  = EEG.srate;
times  = EEG.times;        % times vector (ms)

% Convert times to seconds 
times_s = times / 1000;

%% ---------------------------
% Simple difference kernel 
% ----------------------------
hp_kernel_diff = [-1 1];   % difference kernel (sums to zero -> removes DC)
% Do NOT normalize by sum(abs(...)) for this kernel (keeps intended gain)
% Apply with conv (for demo) and also with filtfilt using an equivalent FIR
hp_conv = conv(signal, hp_kernel_diff, 'same');  % will introduce phase-ish artifacts

% To avoid amplitude scaling confusion, we can scale by 1 (no-op)
% But note this kernel amplifies high-frequency noise; short kernel -> crude HP

%% ---------------------------
%  FIR high-pass filter 
% ----------------------------
% Design parameters
hp_cutoff_hz = 1.0;       % cutoff frequency in Hz (change as needed)
filter_order = 128;       % even, larger -> sharper transition (128 is typical)
if mod(filter_order,2)~=0
    filter_order = filter_order + 1;
end

Wn = hp_cutoff_hz / (srate/2);   % normalized cutoff (0..1, where 1 -> Nyquist)
b_fir = fir1(filter_order, Wn, 'high');  % type: high-pass FIR
a_fir = 1;

% Frequency response for inspection
[H,freqz_w] = freqz(b_fir, a_fir, 1024, srate);

% Apply zero-phase filtering to avoid phase distortion
hp_filtfilt = filtfilt(b_fir, a_fir, signal);

%% PLOTS - kernel and frequency responses

figure('Name','Kernels & Frequency Responses','NumberTitle','off');

subplot(3,2,1)
stem(hp_kernel_diff,'filled')
title('Kernel: simple difference [-1 1]')
xlabel('Samples')

subplot(3,2,2)
plot(freqz_w, abs(H))
xlim([0 srate/2])
xlabel('Frequency (Hz)')
ylabel('Gain')
title(sprintf('FIR HP Response (order=%d, fc=%.1f Hz)', filter_order, hp_cutoff_hz))

% Show FIR kernel (time-domain)
subplot(3,2,3)
plot(b_fir, 'LineWidth', 1.2)
title('FIR kernel coefficients (b\_fir)')
xlabel('Tap index')

% Show zoomed frequency response (low freq)
subplot(3,2,4)
plot(freqz_w, abs(H))
xlim([0 10])
xlabel('Frequency (Hz)')
ylabel('Gain')
title('FIR HP Response (0-10 Hz)')

%% PLOTS - time-domain comparisons
subplot(3,2,5)
plot(times_s, signal, 'k'); hold on
plot(times_s, hp_conv / max(abs(hp_conv))*max(abs(signal)), 'b', 'LineWidth', 1.2); % scaled for visual comparison
legend('Raw','Conv [-1 1] (scaled)')
xlabel('Time (s)')
title('Raw EEG vs conv([-1 1])')

subplot(3,2,6)
plot(times_s, signal, 'k'); hold on
plot(times_s, hp_filtfilt, 'r', 'LineWidth', 1.1)
legend('Raw','FIR hp (filtfilt)')
xlabel('Time (s)')
title('Raw EEG vs FIR HP (zero-phase)')

%% FFT parameters FANCY--------------
N = length(signal);
hz = linspace(0, srate/2, floor(N/2)+1);

fft_raw = fft(signal)/N;
fft_conv = fft(hp_conv)/N;
fft_hpff = fft(hp_filtfilt)/N;

figure('Name','Spectra: Raw vs HP-filtered','NumberTitle','off');
plot(hz, abs(fft_raw(1:length(hz))*2), 'k', 'LineWidth', 1.2); hold on
plot(hz, abs(fft_conv(1:length(hz))*2), 'b', 'LineWidth', 1.0);
plot(hz, abs(fft_hpff(1:length(hz))*2), 'r', 'LineWidth', 1.2);
xlim([0 60])
xlabel('Frequency (Hz)')
ylabel('Amplitude')
legend('Raw','Conv [-1 1]','FIR hp (filtfilt)')
title('Power spectra: raw vs high-pass filtered')

