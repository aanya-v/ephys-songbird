function [reshapedMatrix] = spectsZeroPad(spect_uncut)
% Open compiled spectrograms and pad columns with zero columns to the
% maximum window size of the longest syllable

% Find the maximum number of columns
maxColumns = 0;
for i = 1:length(spect_uncut)
    for j = 1:length(spect_uncut{i})
        currentColumns = size(spect_uncut{i}{j}, 2);
        if currentColumns > maxColumns
            maxColumns = currentColumns;
        end
    end
end

% Pad each cell with zeros to have equal number of columns
for i = 1:length(spect_uncut)
    for j = 1:length(spect_uncut{i})
        currentColumns = size(spect_uncut{i}{j}, 2);
        if currentColumns < maxColumns
            % Padding with zeros
            spect_uncut{i}{j} = [spect_uncut{i}{j}, zeros(98, maxColumns - currentColumns)];
        end
    end
end

% Initialize an empty matrix to hold the reshaped data
numCells = sum(cellfun(@(c) numel(c), spect_uncut)); % total number of cells
reshapedMatrix = zeros(numCells, 98 * maxColumns);

% Counter for the row index in reshapedMatrix
rowIndex = 1;

% Loop through each set and each cell
for i = 1:length(spect_uncut)
    for j = 1:length(spect_uncut{i})
        % Reshape the current cell into a row vector
        currentRow = reshape(spect_uncut{i}{j}, 1, []);

        % Assign it to the reshapedMatrix
        reshapedMatrix(rowIndex, :) = currentRow;
        
        % Update the row index
        rowIndex = rowIndex + 1;
    end
end
end