function [cases_struct] = analy_motor_win_SamAdult(allCases,params_struct)
% DATE: October 2023
% AUTHOR: Leila May Pascual
% PURPOSE:  Compiles syllable-neural cases for spike count/isi analysis of premotor
% windows

% Takes arguments allCases and params_struct
% params_struct must contain:
%       time_window

% OUTPUT FIRING RATE AND ENTROPY VECTORS FOR REAL DATA COMPARISONS
cases_struct = struct();

%%%%   PARAMS    %%%%

% Desired motor window for analysis, in milliseconds
time_window = params_struct.time_window;
casenames= allCases;
case_spikes = {};
ff_case = zeros(length(casenames),1);

for i = 1:length(casenames) % loop through each case
    
    % Load spike trains
    load(strcat(casenames{i},'.mat'),'spiketrains')
    case_trials = spiketrains{2,1};
    case_n_trials(i) = length(case_trials);
    n_spikes = [];

    % count the number of spike counts for each trial
    for j = 1:case_n_trials(i) 
        n_spikes(j) = length(case_trials{j});
    end

    % calculate the average isi for each trial
    mean_isi_case = nan(case_n_trials(i),1);
    mean_ifr_case = nan(case_n_trials(i),1);
    for jj = 1:case_n_trials(i)
        spiketimes= case_trials{jj}; % in milliseconds
        if length(spiketimes) > 1
            mean_isi_case(jj) = mean(diff(spiketimes))/1000; % convert to seconds
            mean_ifr_case(jj) = 1/mean_isi_case(jj);
        end
    end

    % calculate across-trials average spike count & spike rate in Hz
    mean_n_spikes(i) = mean(n_spikes);
    mean_hz(i) = mean_n_spikes(i) / time_window;

    % calculate across-trials averae isi and ifrs and ff
    mean_isi(i) = mean(mean_isi_case,'omitnan'); 
    mean_ifr(i) = mean(mean_ifr_case,'omitnan');
    ff_case(i) = var(n_spikes)/ mean(n_spikes);

    % for each unique spike count, calculate probability
    unique_n = unique(n_spikes);
    prob_vec_case = [];
    for k = 1:length(unique_n)
        % calc the occurence of a given spike rate out of n total trials
        prob_vec_case(k)=length(find(n_spikes==unique_n(k)))/case_n_trials(i);
    end
        
    % calculate the entropy of spike rates for the given spike rate
    entropy_case_adult(i)=-sum(prob_vec_case.*log2(prob_vec_case));

end

% Output Struct file
cases_struct.mean_fr(:,1) = mean_hz;
cases_struct.entropy(:,1) = entropy_case_adult;
cases_struct.case_sz(:,1) = case_n_trials;
cases_struct.casenames(:,1) = casenames;
cases_struct.mean_isi(:,1) = mean_isi;
cases_struct.mean_ifr(:,1) = mean_ifr;
cases_struct.fano_factor(:,1) = ff_case;
end