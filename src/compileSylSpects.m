function [spects_all,labels_all,days_all] = compileSylSpects(birdname,syls,days,savedir,spect_params)
%   Original write date: Dec 2023
%   Author: Leila May Pascual

% Initialize feature matrix
spects_all = [];
labels_all = [];
days_all = [];
labels_by_day = {};
spect_uncut = {};
raw_audio = {};
failedSyls ={};

% print spect parameters
spect_params

% Concatenate spectrograms, syllable labels, and days
for d= 1:length(days)
    day = days(d)
    % change to folder containing files
    [~, ~,~,directory,~] = raster_params(birdname,day);
    cd(directory)

    spect_day = [];
    labels_day = [];
    spect_uncut_day=[];
    raw_audio_day = {};
    try
    [spect_day, labels_day,spect_uncut_day, raw_audio_day] = getSpectMatrix(syls, spect_params);
    [spects_all] = [spects_all; spect_day];
    [labels_all] = [labels_all, labels_day];
    [days_all] = [days_all; repmat(day,size(spect_day,1),1)];
    spect_uncut{d} = spect_uncut_day;
    labels_by_day{d} = labels_day;
    raw_audio = [raw_audio;raw_audio_day];

    catch exception
        message = ['failed for: ' num2str(day) ' -syl ' syls];
        failedSyls = [failedSyls; message];
        disp(message)
    end

    cd ..

end

% Replace values at indices 661 and 662 with 66
days_all(days_all == 661 | days_all == 662) = 66;
days(days == 661 | days == 662) = 66;

% %% Create a spectrogram matrix of all trials padding with zero trials that are less than the maximum syllable duration
% spects_all_zero_padded = spectsZeroPad(spect_uncut);

%% Save output to processed data folder
cd(savedir)

%Convert array of days to a string with underscores between them
daysString = strjoin(arrayfun(@(x) num2str(x), days, 'UniformOutput', false), '_');
dateGenerated = datetime();

% Concatenate the parts to create the filename and save
whichWinStr = sprintf('_syldur_%d_to_%d', spect_params.whichWindows(1), spect_params.whichWindows(2));
filename = strcat('spectrograms_',birdname,'_syls_', syls, '_days_', daysString, whichWinStr, '.mat');
save(filename,"birdname","syls","days","spects_all","labels_by_day",...
    "labels_all","days_all","spect_uncut", "raw_audio", "spect_params",...
    "dateGenerated","failedSyls",'-v7.3' )
% save(filename,"spects_all_zero_padded", '-append');
end