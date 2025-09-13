clc
clear
close all

% ===== Load data =====
load('sampleEEGdata.mat');
signal = EEG.data(10,:,5); % 10th channel, 5th trial
t = (0:length(signal)-1) / EEG.srate;

% Gaussian kernel
time = -1:1/EEG.srate:1;
s = 5/(2*pi*30);
kernel = exp((-time.^2)/(2*s^2))/30;
kernel = kernel / sum(kernel);

% Convolution length
conv_length = length(signal) + length(kernel) - 1;

%% ====== Padding Variants (manual, no padarray) ======

% 1. Symmetric padding (zeros appended at end for both)
sig_sym = [signal, zeros(1, conv_length - length(signal))];
ker_sym = [kernel, zeros(1, conv_length - length(kernel))];

% 2. Asymmetric padding (zeros only in front of kernel, signal padded at end)
sig_asym1 = [signal, zeros(1, conv_length - length(signal))];
ker_asym1 = [zeros(1, conv_length - length(kernel)), kernel];

% 3. Asymmetric padding (zeros only in front of signal, kernel padded at end)
sig_asym2 = [zeros(1, conv_length - length(signal)), signal];
ker_asym2 = [kernel, zeros(1, conv_length - length(kernel))];

% 4. Double-asymmetric padding (kernel split with pre+post zeros, signal normal padded at end)
pad_left = floor((conv_length - length(kernel))/2);
pad_right = ceil((conv_length - length(kernel))/2);
ker_asym3 = [zeros(1, pad_left), kernel, zeros(1, pad_right)];
sig_asym3 = [signal, zeros(1, conv_length - length(signal))];

%% ====== Convolution in Frequency Domain ======
conv_sym   = ifft(fft(sig_sym)   .* fft(ker_sym));
conv_asym1 = ifft(fft(sig_asym1) .* fft(ker_asym1));
conv_asym2 = ifft(fft(sig_asym2) .* fft(ker_asym2));
conv_asym3 = ifft(fft(sig_asym3) .* fft(ker_asym3));

%% ====== Plot Padded Signals ======
figure;
subplot(4,2,1);
plot(sig_sym); title('Signal - Symmetric Pad (post zeros)');
subplot(4,2,2);
plot(ker_sym); title('Kernel - Symmetric Pad (post zeros)');

subplot(4,2,3);
plot(sig_asym1); title('Signal - Post zeros');
subplot(4,2,4);
plot(ker_asym1); title('Kernel - Pre zeros');

subplot(4,2,5);
plot(sig_asym2); title('Signal - Pre zeros');
subplot(4,2,6);
plot(ker_asym2); title('Kernel - Post zeros');

subplot(4,2,7);
plot(sig_asym3); title('Signal - Post zeros');
subplot(4,2,8);
plot(ker_asym3); title('Kernel - Split zeros (pre+post)');

sgtitle('Different Zero-Padding Schemes (Manual)');

%% ====== Plot Convolution Results ======
figure;
plot(real(conv_sym), 'LineWidth',1.5); hold on;
plot(real(conv_asym1), 'LineWidth',1.2);
plot(real(conv_asym2), 'LineWidth',1.2);
plot(real(conv_asym3), 'LineWidth',1.2);
legend('Symmetric (post-post)', 'Asym1 (sig-post, ker-pre)', ...
       'Asym2 (sig-pre, ker-post)', 'Asym3 (sig-post, ker-split)');
xlabel('Samples'); ylabel('Amplitude');
title('Convolution Results under Different Zero-Padding');

%% ====== Frequency Domain View ======
N_fft = conv_length;
frequencies = (0:N_fft-1) * (EEG.srate/N_fft);

figure;
plot(frequencies, abs(fft(conv_sym)), 'LineWidth',1.5); hold on;
plot(frequencies, abs(fft(conv_asym1)));
plot(frequencies, abs(fft(conv_asym2)));
plot(frequencies, abs(fft(conv_asym3)));
xlim([0 EEG.srate/2])
xlabel('Frequency (Hz)'); ylabel('Magnitude');
legend('Symmetric', 'Asym1', 'Asym2', 'Asym3');
title('FFT of Convolution Results (Effect of Padding)');
