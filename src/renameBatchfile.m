function renameBatchfile(inputFileName,outputFileName)
% takes a .txt file containing a list of filenames with .neuralnot.CHx.mat
% extension and removes the extension to the intan recfile name and creates
% a UTF-8 file

% Open the input file for reading
inputFile = fopen(inputFileName, 'r');
% Open the output file for writing
outputFile = fopen(outputFileName, 'w', 'n', 'UTF-8'); % Write in UTF format

% Check if the files are opened successfully
if inputFile == -1 || outputFile == -1
    error('Error in opening files');
end

% Read each line from the input file, process it, and write to the output file
while ~feof(inputFile)
    line = fgets(inputFile);  % Read a line
    if line == -1
        break;  % Exit the loop if end of file is reached
    end

    % Remove '.neuralnot_CH9.mat' and newline characters from the line
    modifiedLine = strrep(line, '.neuralnot_CH9.mat', '');
    modifiedLine = strrep(modifiedLine, newline, '');
    % Write the modified line to the output file
    fprintf(outputFile, '%s\n', modifiedLine);

end

% Close the opened files
fclose(inputFile);
fclose(outputFile);

disp(['Processed filenames written to ', outputFileName]);

end