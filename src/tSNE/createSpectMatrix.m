function [reshapedMatrix] = createSpectMatrix(spect_uncut,selected_bins)
% Open compiled spectrograms and pad columns with zero columns to the
% desired window size 'maxColumns'

n_time_bins = diff(selected_bins)+1;
n_trials = length(spect_uncut);
freq_binlength = length(spect_uncut{1});

% Pad each cell with zeros to have equal number of columns
for i = 1:n_trials
    currentColumns = size(spect_uncut{i}, 2);
    if currentColumns < n_time_bins
        % Padding with zeros
        spect_uncut{i} = [spect_uncut{i}, zeros(freq_binlength, n_time_bins - currentColumns)];
    end

end

% Initialize an empty matrix to hold the reshaped data
reshapedMatrix = zeros(n_trials, freq_binlength * n_time_bins);

% Counter for the row index in reshapedMatrix
rowIndex = 1;

% Loop through each set and each cell
for i = 1:n_trials
    curr_spect = spect_uncut{i};
    % Select time bins for current spectrogram
    spect_select = curr_spect(:,selected_bins(1):selected_bins(2));
    % Reshape the current cell into a row vector
    currentRow = reshape(spect_select, 1, []);

    % Assign it to the reshapedMatrix
    reshapedMatrix(rowIndex, :) = currentRow;
    
    % Update the row index
    rowIndex = rowIndex + 1;
end
end