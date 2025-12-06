%% Convolution Theorem Visualization
% Convolution in time domain = Multiplication in frequency domain

clear; close all; clc;

%% Define time axis
Fs = 1000;            % sampling frequency (Hz)
T  = 1/Fs;            % sampling period
L  = 512;             % number of samples
t  = (0:L-1)*T;       % time vector (0 to ~0.5s)

%% Define two signals
x = double(t>=0.05 & t<=0.15);   % square pulse
h = exp(-50*t);                  % exponential decay

%% --- TIME DOMAIN CONVOLUTION (integration formula) ---
% MATLAB's conv is equivalent to integral of shifted product
y_time = conv(x,h)*T;            % numerical approx of continuous integral
ty = (0:length(y_time)-1)*T;     % new time axis

%% --- FREQUENCY DOMAIN (FFT method) ---
N = length(x) + length(h) - 1;   % proper length for linear convolution
X = fft(x,N);
H = fft(h,N);
Y = X .* H;                      % multiplication in freq domain
y_freq = ifft(Y);                % back to time domain

%% --- Plotting ---
figure;

subplot(3,2,1);
plot(t,x,'b','LineWidth',1.5);
title('Signal x(t)'); xlabel('Time (s)'); ylabel('x(t)');

subplot(3,2,2);
plot(t,h,'r','LineWidth',1.5);
title('Kernel h(t)'); xlabel('Time (s)'); ylabel('h(t)');

subplot(3,2,3);
plot(ty,y_time,'k','LineWidth',1.5);
title('Time-domain convolution'); xlabel('Time (s)'); ylabel('y(t)');

subplot(3,2,4);
plot(ty,y_freq,'g--','LineWidth',1.5);
title('Frequency-domain convolution'); xlabel('Time (s)'); ylabel('y(t)');

subplot(3,2,5:6);
plot(ty,y_time,'k','LineWidth',1.5); hold on;
plot(ty,y_freq,'g--','LineWidth',2);
legend('Time domain conv','Freq domain conv (IFFT)')
title('Comparison: Both results match'); xlabel('Time (s)');
