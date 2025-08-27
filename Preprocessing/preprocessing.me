%% -------------------------
% Load EEG dataset
% -------------------------

eeglab %% launch toolbox
EEG = pop_loadset('eeglab_data.set', fullfile(fileparts(which('eeglab')), 'sample_data'));
% a sample dataset provided with EEGLAB
data = EEG.data;       % channels x samples
chan_labels = {EEG.chanlocs.labels};  %extract chennel labels 

% -------------------------
% Step 1: Filtering: Notch filter at 60 Hz
% -------------------------

%eeg electrodes being very sensitive can pick up even the electrical power
%supply. A notch filter removes a veery narrow band of frequencies.
wo = 60/(EEG.srate/2); bo = wo/35; % w0-> normalized notch freq (60 Hz, normalising using Nyquist formala- sampling rate/2), b0->bandwidth of the notch filter
[bn,an] = iirnotch(wo, bo); % Design the filter
data_filt = filtfilt(bn,an,data')'; % Apply filter, filtfilt applies the filter forward and backward. Therefore no phase distortion

% -------------------------
% Step 2: Plot PSD before & after filtering
% -------------------------
%PSD-Power Spectral Desity 
figure;
[~,~,~,~,~] = spectopo(data, size(data,2), EEG.srate); %Plot PSD of raw data
figure;
[~,~,~,~,~] = spectopo(data_filt, size(data_filt,2), EEG.srate); %Plot PSD of filtered data

%% -------------------------
% Step 3: Referencing: Common Average Reference
% -------------------------
% Voltage values recorded is a relative measure between the active electrode and reference electrode,
% so when rereferencing to the average of all the electrodes which are also relative to
% the reference electrode the reference arbitary electrode cancels outand
% each electrode value i selatove to the average of all electrodes 
%  
avg_ref = mean(data_filt, 1); %take the mean of the filtered data across channels(1 indicates across rows)
data_car = data_filt - avg_ref;% rereferencing

%% -------------------------
% Step 4: Epoching around "rt"
% -------------------------
event_latencies = [EEG.event.latency]; %  Extract latencies of event
event_types     = {EEG.event.type}; % EEG.event has types- 'rt', 'square'
target_idx = find(strcmp(event_types,'rt')); % FInd indices of events with rt, these are the events we want to epoch
epoch_window = round([-0.2 0.8]*EEG.srate);  % epoch : -200 to +800 ms around an event, multiplying by EEG srate converts seconds to samples, and then round it to the nearest integer
epoch_len = diff(epoch_window)+1; % computes the number of samples in each epoch

% Preallocate
% 3D matrix : channels * time * trials
epochs = nan(size(data_car,1), epoch_len, length(target_idx)); % intialisation
for i = 1:length(target_idx) % the identified event types = rt
    center = round(event_latencies(target_idx(i))); % get the time at which the specific event happened
    idx = center+epoch_window(1):center+epoch_window(2); % adds the relative window to the event time, you will get an epoch window
    if idx(1)>0 && idx(end)<=size(data_car,2) % make sure the requested window fits inside the data
        epochs(:,:,i) = data_car(:,idx); % extract the epochs, and stores it fo the ith trial
    end
end

%% -------------------------
% Step 5: Baseline correction
% -------------------------
baseline_idx = 1:round(0.2*EEG.srate); % In the epochs , the first 200 ms (pre-stimulus) samples are baseline 
baseline = mean(epochs(:,baseline_idx,:),2); % Compute mean baseline per channel, per trial
epochs_bc = epochs - baseline; % Subtract baseline from each trial
% epoch_bc : channels*time*trials

%% -------------------------
% Step 6: Improved Artifact Rejection
% -------------------------
% (a) Reject bad trials
% Computes the variance of each channel within each trial, averages that across channels per trial, 
% and then rejects trials whose average variance is more than 3 standard deviations 
% above the global variance mean (across trials).
trial_var = squeeze(var(epochs_bc,0,2)); % epoch_bc-- baseline corrected epochs (compute variance along dimension 2 ie. across time points)
%squeeze removes any singleton dimension you dont want to work with
trial_reject = (mean(trial_var,1) > mean(mean(trial_var))+3*std(mean(trial_var)));
%  mean(trial_var,1)--average variance across channels per trial gives a row vector, one value
% for each trial.
% (mean(trial_var))--global average across all channels and trials mean
% std(mean(trial_var))--standard deviation of trial variances (across trials).
% reject trials where trial variance > global mean + 3*std
good_trials = ~trial_reject; % keep trials that arent rejected
epochs_clean = epochs_bc(:,:,good_trials); % keep only the good trials; 3rd dim

fprintf('Rejected %d/%d trials\n', sum(trial_reject), length(trial_reject));

% (b) Reject bad channels (flat or high variance across trials)
chan_var = squeeze(var(epochs_clean,0,[2 3])); % compute variance across time[dim2] and trials[dim3]
bad_chans = (chan_var < 1e-6) | (chan_var > mean(chan_var)+3*std(chan_var)); 
%if a chnas variance is near zero then its abad channel
% if a channels variance is grreater teh mean+3*std then its dominated by
% noise and artifacts
% channel is bad if either of the condition is true
good_chans = find(~bad_chans);
epochs_clean = epochs_clean(good_chans,:,:); %keep only the good channels; first dim

fprintf('Rejected %d/%d channels\n', sum(bad_chans), length(chan_labels));

chan_labels_clean = chan_labels(good_chans);
