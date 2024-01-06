%% This function checks for overlaps in songbouts excised from the original 
%% 1-minute intan recording files. 

% Author: Leila May Pascual
% Date: Dec 2023

% Output will be a figure that shows which songbout segments overlap

function [overlaps] = findOverlaps(filelist)

if nargin == 1
%% Find duplications in filenames
% input filelist is a .txt file of filenames
inputList = readlines(filelist,'EmptyLineRule','skip');

% files that will be checked for duplication
recmat = cellstr(inputList);

% Initialize a logical array to track duplicates
isDuplicate = false(size(inputList));

% Iterate through each element in the list
for i = 1:length(inputList)
    % Check if the current element is a duplicate
    isDuplicate(i) = any(strcmp(inputList(i), inputList([1:i-1, i+1:end])));
end

% Check if there are any duplicates
    hasDuplicates = any(isDuplicate);
    
    % Display the result
    if hasDuplicates
        disp('The list contains duplicates:');
        disp(sum(isDuplicate))
    else
        disp('The list does not contain duplicates.');
    end

else

% generate list of songbout filenames in the directory
files = dir('*_songbout*.mat'); files = struct2table(files);
for fn = 1:size(files,1)
    filename = char(files.name(fn));
    isrec(fn,1) = length(filename) < 39 && length(filename) > 35;
end
recmat = files.name(isrec);

end

%% LOAD SONGBOUT AND SYLLABLE ONSET AND OFFSET TIMES ACROSS ALL FILES
% for saving syllable and songobut overlaps
savename = 'syllable_overlaps.mat';
overlaps = struct();


% initialize variables
org_fn = cell(length(recmat),1); % store names of each unique recording file looped through
songbout_file_pre = [];
new_ons_pre = [];
new_offs_pre =[]; 
onsets_all = [];
offsets_all = [];
songbout_ons = zeros(length(recmat),1);
songbout_offs = zeros(length(recmat),1);
songbout_ons_intan = zeros(length(recmat),1);
songbout_offs_intan = zeros(length(recmat),1);


% extract songbout and syllable times from songbout files
for i = 1:length(recmat)
    % open songbout recording file
    filename =recmat{i};
    org_fn{i} = filename(1:24); %remove 'songbout*'
    load(filename,'song_onset', 'song_offset','t_board_adc');

    % Store all songbout on/off times
    songbout_ons(i) = song_onset;
    songbout_offs(i) = song_offset;

    % Calculate intan timestamps of songbout on/off times
    songbout_ons_intan(i) = t_board_adc(1);
    songbout_offs_intan(i) = t_board_adc(end);

    fn_notmat = [recmat{i} '.not.mat'];
    if exist(fn_notmat)
        load(fn_notmat,"onsets","offsets");
    
        % convert syllable onsets/offsets to seconds
        onsets_s = onsets./1000;
        offsets_s = offsets./1000;

        % save intan timestamps of onset/offsets to .not.mat
        onsets_intan = (onsets_s + t_board_adc(1))*1000;
        offsets_intan = (offsets_s + t_board_adc(1))*1000;
        save(fn_notmat,'onsets_intan','offsets_intan','-append')

        % calculate and concatenate the intan timestamp of syl onset and offset times
        songbout_file_pre = [songbout_file_pre; repmat(string(filename),length(onsets_s),1)];
        new_ons_pre = [new_ons_pre; onsets_s + t_board_adc(1)];
        new_offs_pre = [new_offs_pre; offsets_s + t_board_adc(1)];

        % concatenate all onset/offset times in relation to songbout file
        onsets_all = [onsets_all; onsets];
        offsets_all = [offsets_all; offsets];

    end
end


%% FIND OVERLAPS IN SYLLABLES


% Sort syllable trials by onset times
[ons_sorted,indices] = sort(new_ons_pre);
new_ons = new_ons_pre(indices);
new_offs = new_offs_pre(indices);
songbout_file = songbout_file_pre(indices);
onsets_all = onsets_all(indices);
offsets_all = offsets_all(indices);

% Check for overlap in time intervals
overlapSyls = new_ons(2:end) < new_offs(1:end-1);
overlaps.syllables_times(:,1) = new_ons(overlapSyls);
overlaps.syllables_times(:,2) = new_offs(overlapSyls);
overlaps.syllables_times(:,3) = new_ons(find(overlapSyls)+1);
overlaps.syllables_times(:,4) = new_offs(find(overlapSyls)+1);
overlaps.syllables_filename(:,1) = songbout_file(overlapSyls);
overlaps.syllables_filename(:,2) = songbout_file(find(overlapSyls)+1);

overlaps.syllables_ons_offs_by_songboutfile(:,1) = songbout_file(overlapSyls);
overlaps.syllables_ons_offs_by_songboutfile(:,2) = onsets_all(overlapSyls);
overlaps.syllables_ons_offs_by_songboutfile(:,3) = offsets_all(overlapSyls);
overlaps.syllables_ons_offs_by_songboutfile(:,4) = songbout_file(find(overlapSyls)+1);
overlaps.syllables_ons_offs_by_songboutfile(:,5) = onsets_all(find(overlapSyls)+1);
overlaps.syllables_ons_offs_by_songboutfile(:,6) = offsets_all(find(overlapSyls)+1);
hasOverlap = sum(new_ons(2:end) < new_offs(1:end-1));

if hasOverlap
    disp('Segments contain syllable overlaps.');
    fprintf('overlaps = %d ', hasOverlap);
else
    disp('Segments do not contain syllable overlaps.');
end

% VISUALIZE SYLLABLE OVERLAPS
% Plot segments
figure;
hold on;

% Plot segments as vertical lines
for i = 1:length(new_ons)
    plot([new_ons(i), new_offs(i)], [i, i], 'LineWidth', 2,'Color','black');
end

% Highlight overlaps in red
if hasOverlap
    for i = 1:length(new_ons)-1
        if new_ons(i+1) < new_offs(i)
            overlapStart = new_ons(i+1);
            overlapEnd = min(new_offs(i), new_offs(i+1));
            plot([overlapStart, overlapEnd], [i+10, i+10], 'r', 'LineWidth', 10);
        end
%         legend(string(overlaps.syllables_times),"Location","southoutside","Interpreter",'latex')
    end
    title('Segments with Syllable Overlaps (Red)');
    subtitle(['Number of Overlaps = ', num2str(hasOverlap)])
else
    title('Segments without Syllable Overlaps');
end

xlabel('Time froms start of Intan recording (s)');
ylabel('Vocal unit (syllable)');
ylim([0, length(new_ons)+1]);
grid on;
hold off;
            

%% FIND OVERLAPS in SONGBOUT FILES
% uses timestamps where t0 is the start of songbout file

% Get the number of files
numFiles = numel(org_fn);

% Initialize a logical array to store overlap information
overlapMatrix = false(numFiles);

% Iterate through each pair of files
for i = 1:numFiles
    for j = i+1:numFiles
        % Check if the filenames are the same
        if strcmp(org_fn{i}, org_fn{j})
            % Check for overlap in time intervals
            overlap = (songbout_ons(i) < songbout_offs(j)) && (songbout_offs(i) > songbout_ons(j));
            
            % Store the overlap information in the matrix
            overlapMatrix(i, j) = overlap;
            overlapMatrix(j, i) = overlap;
        end
    end
end


% Find rows with overlaps
filesWithOverlaps = any(overlapMatrix, 2);
org_fn_overlaps = org_fn(filesWithOverlaps);
files_overlaps = recmat(filesWithOverlaps);
song_ons_overlaps = songbout_ons(filesWithOverlaps);
song_offs_overlaps = songbout_offs(filesWithOverlaps);

overlaps.songbouts(:,1) = string(files_overlaps);
overlaps.songbouts(:,2) = song_ons_overlaps;
overlaps.songbouts(:,3) = song_offs_overlaps;

% Find overlaps using timestamps from intan

    % sort by sonbout onset
[songbout_ons_intan_sorted, ind] = sort(songbout_ons_intan);
songbout_ons_intan_sorted = songbout_ons_intan(ind);
songbout_offs_intan_sorted = songbout_offs_intan(ind);
recmat_sorted = recmat(ind);
overlapSongbouts = songbout_ons_intan_sorted(2:end) < songbout_offs_intan_sorted(1:end-1);
overlapSongboutsFn = recmat_sorted(overlapSongbouts);
overlaps.songbouts_intan(:,1) = string(overlapSongboutsFn);
overlaps.songbouts_intan(:,2) = songbout_ons_intan_sorted(overlapSongbouts);
overlaps.songbouts_intan(:,3) = songbout_offs_intan_sorted(overlapSongbouts);
hasOverlapSongbouts = sum(overlapSongbouts);

if hasOverlapSongbouts
    disp('Found songbout overlaps.');
    fprintf('Songbout overlaps = %d', hasOverlapSongbouts);
else
    disp('No songbout overlaps found.');
end


%% VISUALIZE OVERLAPS

figure; clf;
% Initialize a figure with a larger height
% figure('Position', [100, 100, 800, 900]);
% Save current warning state
warning('off')
% Define colors for identical filenames
uniqueFilenames = unique(org_fn_overlaps);
numUniqueFilenames = numel(uniqueFilenames);
colorMap = turbo(numUniqueFilenames);

% Create a colormap for filenames
filenameColors = zeros(length(files_overlaps), 3);  % Initialize the color matrix

% Assign colors to filenames
for i = 1:length(org_fn_overlaps)
    colorIndex = find(strcmp(uniqueFilenames,org_fn_overlaps{i}));
    filenameColors(i, :) = colorMap(colorIndex, :);
end

if length(song_ons_overlaps) > 1
    for j = 1:length(song_ons_overlaps)
        plot([song_ons_overlaps(j), song_offs_overlaps(j)], [-j, -j], 'Color',filenameColors(j,:), 'LineWidth', 2);
        hold on;
    end

    ylim([-length(song_ons_overlaps)-1 0])    
end
title('Overlapping songbout segments found')
subtitle('empty plot means no songbout overlaps found')
legend(files_overlaps,'Location', 'eastoutside','Interpreter','latex')


% Save duplications
save(savename,'recmat', 'songbout_ons','songbout_offs','overlaps',...
    'songbout_ons_intan', 'songbout_offs_intan')
end