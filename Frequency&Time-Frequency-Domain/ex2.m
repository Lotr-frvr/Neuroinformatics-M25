%% Load 50 points of EEG data (from one electrode)
% For demo, we'll just simulate EEG as filtered noise
% Replace this line with your real EEG vector (length = 50)
EEGdata = filter(1,[1 -0.95],randn(1,50));  

time = 1:50;

%% Define kernels
kernel_U     = [0 1 2 1 0];
kernel_U     = kernel_U / sum(kernel_U);    % normalize
kernel_decay = [1 0.5 0.25 0.125];
kernel_decay = kernel_decay / sum(kernel_decay);

%% Convolve manually (time-domain convolution)
conv_U     = conv(EEGdata, kernel_U, 'same');
conv_decay = conv(EEGdata, kernel_decay, 'same');

%% Plot
figure;

subplot(3,1,1)
plot(time, EEGdata, 'k', 'LineWidth',1.5)
title('EEG data (50 time points)')
xlabel('Time (samples)'), ylabel('Amplitude')

subplot(3,1,2)
plot(kernel_U, 'bo-','LineWidth',1.5); hold on
plot(kernel_decay, 'ro-','LineWidth',1.5);
legend('Inverted U kernel','Decay kernel')
title('Kernels')

subplot(3,1,3)
plot(time, EEGdata, 'k','LineWidth',1.5); hold on
plot(time, conv_U, 'b','LineWidth',2)
plot(time, conv_decay, 'r','LineWidth',2)
legend('Original EEG','Convolved w/ U','Convolved w/ Decay')
title('EEG convolved with kernels')
xlabel('Time (samples)')
