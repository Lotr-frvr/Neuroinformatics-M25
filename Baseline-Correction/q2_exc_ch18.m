% Select five time points and create topographical maps of power with
% and without baseline normalization at each selected time-frequency
% point. You should have time in columns and with/without baseline
% normalization in rows. Use separate figures for each frequency. The
% color scaling should be the same for all plots over time within a
% frequency, but the color scaling should be different for with versus
% without baseline normalization and should also be different for each
% frequency.


clear; clc;
load sampleEEGdata

%% --- wavelet parameters----------------------------------
min_freq = 2;
max_freq = 40;      % Reduce for speed; change to 128 if needed
frequencies = [2 10 15 20 25 30 40];   % Hz
num_frex = length(frequencies);

time = -1:1/EEG.srate:1;
half_wav = (length(time)-1)/2;
wavelet_cycles = 4;

n_wavelet = length(time);
n_data    = EEG.pnts;
n_conv    = n_wavelet+n_data-1;
n_pow2    = pow2(nextpow2(n_conv));

%% --- Prepare TF matrix: channels × freqs × time -------------------------
tf_raw = zeros(EEG.nbchan, num_frex, EEG.pnts);

for ch = 1:EEG.nbchan
    
    % FFT of channel data (averaged over trials)
    chan_data = squeeze(mean(EEG.data(ch,:,:),3));
    fft_data = fft(chan_data, n_pow2);
    
    for fi = 1:num_frex
        
        f = frequencies(fi);
        s = wavelet_cycles / (2*pi*f);
        
        % Complex Morlet wavelet
        wavelet = exp(2*1i*pi*f.*time) .* exp(-time.^2 ./ (2*s^2));
        fft_wave = fft(wavelet, n_pow2);
        
        % Convolution
        conv_res = ifft(fft_wave .* fft_data);
        conv_res = conv_res(1:n_conv);
        conv_res = conv_res(half_wav+1:end-half_wav);
        
        tf_raw(ch,fi,:) = abs(conv_res).^2;
    end
end


%% --- baseline normalization (dB) ----------------------------------------
baseline_ms = [-500 -200];
[~,b1] = min(abs(EEG.times-baseline_ms(1)));
[~,b2] = min(abs(EEG.times-baseline_ms(2)));

baseline_power = squeeze(mean(tf_raw(:,:,b1:b2),3));   % chan × freq
tf_db = 10*log10( bsxfun(@rdivide, tf_raw, baseline_power) );


%% --- Pick 5 time points -------------------------------------------------
timepoints_ms = [-100 0 200 400 800];
time_idx = zeros(size(timepoints_ms));
for i=1:length(timepoints_ms)
    [~,time_idx(i)] = min(abs(EEG.times - timepoints_ms(i)));
end


%% --- Topographical plots ------------------------------------------------
% Uses the EEGLAB function topoplot
% Make sure EEGLAB is in your MATLAB path

for fi = 1:num_frex
    
    f = frequencies(fi);
    
    % Extract topo data for the 5 time points
    raw_maps = squeeze(tf_raw(:,fi,time_idx));
    db_maps  = squeeze(tf_db(:,fi,time_idx));
    
    % Compute color scales
    raw_clim = [min(raw_maps(:)) max(raw_maps(:))];
    db_clim  = [min(db_maps(:))  max(db_maps(:))];
    
    figure('Name',sprintf('Frequency %.1f Hz',f),'Position',[100 100 1200 600]);
    
    % ----- Row 1: RAW power -----
    for t = 1:length(time_idx)
        subplot(2, length(time_idx), t)
        topoplot(raw_maps(:,t), EEG.chanlocs, 'maplimits', raw_clim);
        title(sprintf('%d ms (RAW)', timepoints_ms(t)))
        colorbar
    end
    
    % ----- Row 2: dB normalized -----
    for t = 1:length(time_idx)
        subplot(2, length(time_idx), length(time_idx)+t)
        topoplot(db_maps(:,t), EEG.chanlocs, 'maplimits', db_clim);
        title(sprintf('%d ms (dB)', timepoints_ms(t)))
        colorbar
    end
    
    sgtitle(sprintf('Topographical Maps at %.2f Hz', f));
end

