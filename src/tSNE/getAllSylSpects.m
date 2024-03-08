function [spectrograms_full_syl_durations,labels_all,days_all] = getAllSylSpects(savedir,output_filename,params)
%   Original write date: Jan 2023
%   Evolved from getSylSpects
%   Author: Leila May Pascual

% Initialize output variables
spectrograms_full_syl_durations = {};
labels_all = [];
days_all = [];
raw_audio = {};
failedSyls ={};

% print spect parameters
params
birdname = params.birdname;
syls = params.syls;
days = params.days_folder_id;

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
    [spect_day, labels_day, raw_audio_day] = getSpects(syls, params);

    % concatenate trials
    spectrograms_full_syl_durations = [spectrograms_full_syl_durations; spect_day];
    labels_all = [labels_all, labels_day];
    days_all = [days_all; repmat(day,size(spect_day,1),1)];
    raw_audio = [raw_audio;raw_audio_day];

    catch exception
        message = ['failed for: ' num2str(day) ' -syl ' syls];
        failedSyls = [failedSyls; message];
        disp(message)
    end

end

% Replace values at indices 661 and 662 with 66
days_all(days_all == 661 | days_all == 662) = 66;
days(days == 661 | days == 662) = 66;
days_all(days_all == 801 | days_all ==802 | days_all ==803) = 80;
days(days ==801 | days ==802 | days==803) = 80;
days_all(days_all == 811 | days_all ==812) = 81;
days(days ==811 | days ==812) = 81;

%% Save output to processed data folder
dateGenerated = char(datetime("today","Format","uuuu-MMM-dd"));
cd(savedir)
save(output_filename,"birdname","syls","days","params",...
    "labels_all","days_all","spectrograms_full_syl_durations", "raw_audio", ...
    "dateGenerated","failedSyls",'-v7.3' )

end