function [num_time_bins] = calcNumWindows(audio_length_ms,overlap_windowsize)

% Parameters
FS = 30000; % Sampling rate in Hz
spect_win_dur = overlap_windowsize(2); % Window size in milliseconds
spect_overlap = overlap_windowsize(1); % Pct overlap

% Convert milliseconds to samples
window_size_samples = FS * spect_win_dur / 1000; % Window size in samples
audio_length_samples = FS * audio_length_ms / 1000; % Audio length in samples

% Calculate step size (window shift per time bin) in samples
overlap_samples = window_size_samples * spect_overlap; % Overlap in samples
step_size_samples = window_size_samples - overlap_samples; % Step size in samples

% Calculate the number of time bins
num_time_bins = ceil((audio_length_samples - window_size_samples) / step_size_samples + 1);

% Display the result
fprintf('Number of time bins: %d\n', num_time_bins);


end