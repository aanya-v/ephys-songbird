function [spects_matrix,labels_all] = getSpectMatrix(syllables,params)
% Import the necessary libraries
import matlab.io.*
import statistics.toolbox.*

%% Specify spectogram parameters
Fs = 30000; %sampling rate in Hz

spect_params = params.overlap_windowsize;
numWindows = params.numWindows;
if nargin ==1    
    spect_params = [0.9 4]; %[percent_overlap window_size_in_milliseconds]
    numWindows = 60; % takes the last 60 windows of the spectrogram
end

%% compile spectrogram matrix of syllables
%  load sound vectors and syllable onset/offset times
soundfiles = dir('*_songbout*.not.mat'); files = struct2table(soundfiles);
for fn = 1:size(files,1)
    filename = char(files.name(fn));
    recname(fn,:) = filename(1:end-8);
    % isrec(fn,1) = length(filename) < 39 && length(filename) > 35;
end
recmat = string(recname);

% Initialize variables
org_fn = cell(length(recmat),1); % store names of each unique recording file looped through
songbout_fn = [];
songbout_ons = zeros(length(recmat),1);
songbout_offs = zeros(length(recmat),1);
songbout_ons_intan = zeros(length(recmat),1);
songbout_offs_intan = zeros(length(recmat),1);
ons_intan_all = [];
offs_intan_all = [];
labels_all = [];
spects_matrix = [];

%% extract songbout and syllable times from songbout files
for i = 1:length(recmat)
    % open songbout recording file
    filename =char(recmat(i));
    underscores = strfind(filename, '_');
    org_fn{i} = [filename(1:underscores(end)-1) '.mat']; %remove 'songbout*'
    load(filename,'song_onset', 'song_offset','t_board_adc','board_adc_data');

    % Store all songbout on/off times
    songbout_ons(i) = song_onset;
    songbout_offs(i) = song_offset;

    % Calculate intan timestamps of songbout on/off times
    songbout_ons_intan(i) = t_board_adc(1);
    songbout_offs_intan(i) = t_board_adc(end);
    
    % store sound vector
    rawsong = board_adc_data;
   

    fn_notmat = [recmat{i} '.not.mat'];
    if exist(fn_notmat)
        load(fn_notmat,"onsets","offsets","labels");
        ons_intan = [];
        offs_intan =[]; 
        select_idx = [];
        labels_select =[];
        
        % choose trials for only desired syllables
        select_idx = ismember(labels, syllables);
        labels_select = labels(select_idx);
        onsets = onsets(select_idx);
        offsets = offsets(select_idx);

        % store all of the labels
        labels_all = [labels_all, labels_select];

        % convert syllable onsets/offsets to seconds
        onsets_s = onsets./1000;
        offsets_s = offsets./1000;

        % calculate and concatenate the intan timestamp of syl onset and offset times
        songbout_fn = [songbout_fn; repmat(string(filename),length(onsets_s),1)];
        ons_intan = onsets_s + t_board_adc(1);
        offs_intan = offsets_s + t_board_adc(1);
        ons_intan_all = [ons_intan_all; ons_intan];
        offs_intan_all = [offs_intan_all; offs_intan];

        % get sound vector for each desired syllable
        for x=1:length(labels_select)
            F1 = []; T=[]; P1=[]; 

            on_id = find(abs((t_board_adc-ons_intan(x)))==min(abs(t_board_adc-ons_intan(x)))); %which time in RAWSONG is closest to syllable onset?
            off_id = find(abs((t_board_adc-offs_intan(x)))==min(abs(t_board_adc-offs_intan(x)))); %...to syllable offset?

            syl_wav = rawsong(on_id:off_id); %raw waveform for current syllable
            [S1] = spect_from_waveform(syl_wav,Fs,0,spect_params); %get syllable spectrogram
            S{x,i} = S1;
            % cut the spectrograms to equal duration lengths
            S_cut = S1(:,end-numWindows:end);
            % reshape multidim spec to 1 row for each trial            
            S_collapsed = reshape(S_cut,1,[]);
            Slog = log(1+abs(S_collapsed)); % use log-scale for visualization
            
            % concatenate all spect trials
            spects_matrix = [spects_matrix; Slog]; % use log-scale for visualization
        end
    end
end
end