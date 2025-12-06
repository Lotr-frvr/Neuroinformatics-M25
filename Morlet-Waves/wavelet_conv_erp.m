
%% 
% Average the result of convolution over all trials and plot an ERP
% corresponding to each wavelet frequency. Each frequency should be in
% its own subplot.

% Average across trials (ERP per frequency)`
% 
% What it does:
% 1. Builds a smaller set of wavelets (e.g., 10, 2–80 Hz, log-spaced).
% 2. Concatenates trials into one long vector and computes a single FFT-based convolution per wavelet.
% 3. Reshapes the result to (time × trials), then averages to obtain ERPs per frequency.


%% Q1: Generate wavelets 
clc; close all;

load sampleEEGdata

min_freq = 2;
max_freq = 80;
num_frex = 10;

frex = logspace(log10(min_freq), log10(max_freq), num_frex); % just 80-2 diving into 30 steps
time = -1:1/EEG.srate:1;

% width (same as reference code)
s    = logspace(log10(3),log10(10),num_frex) ./ (2*pi*frex);

wavelets = cell(num_frex,1);

for fi = 1:num_frex
    wavelets{fi} = exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(s(fi)^2)));
end

%% Q2: Convolution with all trials 


chan2use = 'fcz';
data = squeeze(EEG.data(strcmpi(chan2use,{EEG.chanlocs.labels}),:,:));
data = reshape(data,1,[]);   % all trials concatenated

n_data = EEG.pnts * EEG.trials;

eegconv_all = cell(length(frex),1);

for fi = 1:length(frex)
    
    wavelet = wavelets{fi};
    n_wavelet = length(wavelet);
    n_conv = n_wavelet + n_data - 1;

    % FFTs
    fft_wavelet = fft(wavelet, n_conv);
    fft_data = fft(data, n_conv);

    % convolution
    eegconv = ifft(fft_wavelet .* fft_data);
    
    halfw = ceil(n_wavelet/2);
    eegconv = eegconv(halfw:end-halfw+1);

    % reshape back to time × trials
    eegconv_all{fi} = reshape(eegconv, EEG.pnts, EEG.trials);
end

%% Q3: ERP per wavelet frequency


figure
num_frex = length(frex);

for fi = 1:num_frex
    
    % Extract convolution result for this frequency
    data_fi = eegconv_all{fi};
    
    % ERP = mean real part over trials
    erp = mean(real(data_fi), 2);
    
    subplot(num_frex,1,fi)
    plot(EEG.times, erp, 'k')
    ylabel([num2str(round(frex(fi))) ' Hz'])
    
    if fi < num_frex
        set(gca,'xtick',[])
    else
        xlabel('Time (ms)')
    end
end

sgtitle('ERP per Wavelet Frequency')
