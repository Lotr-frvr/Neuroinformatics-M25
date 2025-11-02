% Morlet Wavelet Family (2 Hz to 30 Hz in 5 steps)
clear; clc; close all;

fs = 200;                % sampling rate (Hz) – adjust as needed
t = -1:1/fs:1;           % time vector (2 sec window)

% Frequencies for the wavelet family
freqs = linspace(2,30,5);

% Gaussian width (standard deviation)
sigma = 0.2;  % wider = smoother envelope

figure;
for i = 1:length(freqs)
    f = freqs(i);

    % Morlet wavelet = Gaussian * complex sine
    wavelet = exp(-t.^2/(2*sigma^2)) .* exp(1i*2*pi*f*t);

    % Plot real part
    subplot(length(freqs),1,i);
    plot(t, real(wavelet), 'b', 'LineWidth', 1.5); hold on;
    plot(t, imag(wavelet), 'r--', 'LineWidth', 1);
    title(['Morlet Wavelet at ' num2str(f) ' Hz']);
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;
end

legend('Real part','Imag part');
sgtitle('Family of Morlet Wavelets (2–30 Hz)');

%%
% Select one electrode from the scalp EEG dataset and convolve each
% wavelet with EEG data from all trials from that electrode. Apply the
% Matlab function real to the convolution result, as in
% convol_result=real(convol_result). This will return the EEG data
% bandpass filtered at the peak frequency of the wavelet. You learn more
% about why this is in the next chapter.
clc 
close all
clear all
% Load EEG data
load sampleEEGdata.mat
% EEG structure is usually called EEG with fields: EEG.data (chan x time x trials), EEG.srate, EEG.times, etc.

fs = EEG.srate;          % sampling rate
times = EEG.times / 1000; % convert ms to sec if needed

% Select one electrode (e.g., Cz, or index 47 in this dataset)
chanIdx = find(strcmpi({EEG.chanlocs.labels}, 'Cz')); % or set manually
signal = squeeze(EEG.data(chanIdx,:,:));  % time x trials

% Define Morlet wavelet family
freqs = linspace(2,30,5);   % 5 wavelets, 2–30 Hz
t = -2:1/fs:2;              % 4-second window
sigma = 0.3;                % Gaussian width

% Loop over wavelets and convolve
convol_results = cell(length(freqs),1);

for fi = 1:length(freqs)
    f = freqs(fi);
    
    % Morlet wavelet
    wavelet = exp(-t.^2/(2*sigma^2)) .* exp(-1i*2*pi*f*t);
    
    % Initialize storage
    nTrials = size(signal,2);
    conv_data = zeros(size(signal));
    
    for tr = 1:nTrials
        % Convolution (same length as data)
        convsig = conv(signal(:,tr), wavelet, 'same');
        % Take the real part (band-pass filtered EEG at f Hz)
        conv_data(:,tr) = real(convsig);
    end
    
    convol_results{fi} = conv_data;
    
    fprintf('Done with %g Hz\n', f);
end

%% Plot an example trial
trialNum = 1;
figure;
for fi = 1:length(freqs)
    subplot(length(freqs),1,fi);
    plot(times, convol_results{fi}(:,trialNum));
    title([num2str(freqs(fi)) ' Hz filtered EEG (real part)']);
    xlabel('Time (s)'); ylabel('Amplitude');
end
sgtitle(['Electrode ' EEG.chanlocs(chanIdx).labels ' - Trial ' num2str(trialNum)]);

erps = cell(length(freqs),1);

for fi = 1:length(freqs)
    conv_data = convol_results{fi};   % time × trials
    erps{fi} = mean(conv_data, 2);    % average over trials
end

 % Plot ERP per frequency
figure;
for fi = 1:length(freqs)
    subplot(length(freqs),1,fi);
    plot(times, erps{fi}, 'k', 'LineWidth', 1.5);
    title([num2str(freqs(fi)) ' Hz ERP (real part)']);
    xlabel('Time (s)'); ylabel('Amplitude');
    grid on;
end
sgtitle(['ERP at Electrode ' EEG.chanlocs(chanIdx).labels]);


%% Broadband ERP (no convolution)
% Average raw signal across trials
broadbandERP = mean(signal,2);   % time x 1

%Plot broadband + wavelet-convolved ERPs
figure;

% 1st subplot: broadband ERP
subplot(length(freqs)+1,1,1);
plot(times, broadbandERP, 'k', 'LineWidth', 1.5);
title(['Broadband ERP at ' EEG.chanlocs(chanIdx).labels]);
xlabel('Time (s)'); ylabel('Amplitude');
grid on;

% Remaining subplots: wavelet ERPs
for fi = 1:length(freqs)
    subplot(length(freqs)+1,1,fi+1);
    plot(times, erps{fi}, 'b', 'LineWidth', 1.5);
    title([num2str(freqs(fi)) ' Hz ERP (real part)']);
    xlabel('Time (s)'); ylabel('Amplitude');
    grid on;
end

sgtitle('Broadband vs. Wavelet-convolved ERPs');

