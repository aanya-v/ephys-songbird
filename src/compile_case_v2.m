%%% DESCRIPTION %%%

%% UPDATES BY LEILA %%
%  Last update: Oct 14 2023
%  This script adapts get_motif_spikes_from_list into creating data
%  structures that is compatible with Sam's adult RA data structures

% FEB 21 2024 - This v2 changes the postmotor cushion to 0.6 sec. 
% It also attempts to allow a sequence of syllables but aligns only to the
% first syllable in the sequence.


% This function takes a file list as input of recording files that 
% passes SAUCY's single unit test and collates neural data during 
% a specified song syllable or sequence.

% update from get_motif_spikes2021. This now adds a postmotor "cushion" of time after the syllable/sequence. 
%This is hard-coded near the top of the main script, and the value is saved in the output as a variable.

% Unlike previous versions (pre-2021) of this function, this does NOT align times to a
% particular point in the motif.

% The final output will be a file containing a data structure called
% "clean_trials" that will include spike times during the specified
% sequence, as well as information about previous and following syllables
% and syllable quantification (pitch, etc) if that is available.

%%%% BEFORE RUNNING THIS SCRIPT:

% 1) Sort spikes in SAUCY

% 2) Label song in evsonganaly

% 3) Optional: Run compare_batchmode_passing_files if you only would like
% to use a subset of files that have simultaneous SUs

% 4) Optional: Run latest version of mark_artifact_periods if necessary (designed to work with v5)

% 5) Optional: Run quantify_pitch_intan_leila

% 6) Optional: Run QuantPitchContours_2020

% Inputs: (1) unit ID without channel number (2) syllable or sequence, (3) neural channel, (4) pre-onset time in seconds,
% (5) compare_batchmode output or batchmode .txt
% output

%% Example syntax: get_motif_spikes_2022('bl136yw14_200126_dv3410','abc',2,.02,'files_with_SUs_ch7_ch14_2020-05-27.mat')

%% Output: A file named
% "[unitID]ch[#]_spiketimes_for_seq_[sequence]_[date]".
% containing information about the inputs to the function, as well as two
% arrays called "clean_trials" and "artifact_trials", based on the
% prescence of artifact defined by mark_artifact_periods. Each array column
% is as follows:
% (1) Originating cbin/rhd
% (2) Syllable onset and offset times within original cbin for this trial
% (not including pre-time)
% (3) syllables before specified sequence (as dictated by get_N_back_syls)
% (4) syllables after sequence (as dictated by get_N_forward_syls)
% (5) Concatenated pre-onset, sequence, and postmotor cushion spiketimes
% (6) Pitches for syllables in specified sequence for that iteration (cell array)
% (7) Entropies (cell array)
% (8) Amplitudes (cell array)
% (9) Times of song quantification (cell array)
% (10) Pitch contours (cell array)
% (11) 


% Other functions called: evsoundin

function [neuralcase] = compile_case(birdID,neuronID,syl_or_seq,neural_channel,pre_onset_time_s,filelist)

%%% MAIN SCRIPT %%%

 tic % start timer

%% SET CRITERIA FOR RETRIEVING TRIALS %%
max_artifact_samples = 20; % Maximum allowable samples containing artifact, as labeled in
% mark_artifact_periods. Not necessarily continuous. One sample = about 0.032 or 0.030 milliseconds (cbin and rhd,
% respectively)

max_intersyl_s = 1; % in s. of any inter-syllable intrevals within the sequence are greater than this the trial will be skipped.
postmotor_cushion = .6; %extra time added to end of sequence in seconds
get_N_back_syls = 6; % tracks syllables and syllable times prior to SYLLABLE OR SEQUENCE entered
get_N_forward_syls = 6; % same as above but following
syl_or_seq = char(syl_or_seq);
max_trials = 1500;

 %% CREATE VARIABLES FOR OUTPUT FILE %%
 case_ID = strcat(birdID,'_',syl_or_seq,'_',neuronID);
 matname = strcat(case_ID, '.mat');
 birdname = birdID;
 syllablename = syl_or_seq;
 neuronname = neuronID;

 recfilename = strings();
 syl_ons_offs = [];
 syl_dur = [];
 syls_prev = strings(max_trials,get_N_back_syls);
 syls_next = strings(max_trials,get_N_forward_syls);
 syl_in_seq = strings(max_trials,1);
 pre_syl_gap =[];
 post_syl_gap = [];
 spiketrains = {};
 acoustics = {};
 date_compiled = [];

neuralcase = struct();
neuralcase.case_ID = case_ID;
neuralcase.birdname = birdID;
neuralcase.syllablename = syllablename;
neuralcase.neuronname = neuronname;

motif_spikes = cell(1,11);
artifact_trials = {};
clean_trials = {};
clean_trials_table = table();

 
%% Load list of recording files to take trials from

if contains(filelist,'batchmoderesults') %% THIS ISN"T WORKING FIX
    batchmode_file = filelist;
    batchmode_data = readtable(batchmode_file,"MissingRule","omitrow",'VariableNamingRule','preserve');
    files = batchmode_data.Var1;
    files = sort(files); % sort files chronologically (alphabetically)

elseif contains(filelist,'txt') % text file should consist only of filenames of passing files
    passedsaucy = readlines(filelist);
    notmats = dir('*not.mat'); notmats = {notmats.name}; 
    % extract only the recname (get rid of other extensions)
      passedsaucy = char(passedsaucy); notmats = char(notmats);
      idx = strfind(passedsaucy(1,:),'mat');
      passedsaucy = passedsaucy(:,1:idx(1)+2); notmats = notmats(:,1:idx(1)+2);
      passedsaucy = string(passedsaucy); notmats = string(notmats);
    % find notmats files that are in the passed list
      files = passedsaucy(ismember(passedsaucy,notmats));
      files = sort(files);

elseif contains(filelist,'files_with_SU')
    load(filelist)
    files = sort(shared_files);

elseif contains(filelist, 'mat')
    if contains(filelist, 'clust3')
        load(filelist,'passed_clust3')
        files = passed_clust3;
    elseif contains(filelist, 'clust2')
        load(filelist,'passed_clust2')
        files = passed_clust2;
    end
else
    error('File list should either be a batchmode ouput .txt or compare_batchmode .mat file. Check that the filename was entered correctly,')
end


if contains(files{1},'rhd') % check file type
     filetype = 'rhd';
elseif contains(files{1},'mat')
    filetype = 'mat';
elseif contains(files{1},'cbin')
    filetype = 'cbin';
else
    error('Something went wrong; first file on the list is not a .cbin or .rhd or .mat')
end


%% loop through all files in the list to retrieve syllable trials
   ct = 0;
   for i = 1:size(files,1)

       current_fname = files{i,1};
       if contains(current_fname, 'neuralnot')
           current_fname = current_fname(1:end-19);
       end

       clearvars weighted_avg_labelvec pitchContours
       % load first set of syllable (notmat) and sorted spike (neuralnot) files
       load([current_fname(1:end-4),'.mat.not.mat'],'labels','onsets','offsets','weighted_avg_labelvec',...
           'wiener_entropy_labelvec','amp_at_pitchquant_labelvec','t_assay_labelvec','pitchContours',...
           'spect_entropy_labelvec') % song labels and syllable times and acoustics
       load([current_fname,'.neuralnot_CH',num2str(neural_channel),'.mat']) % spike sorting data
       
       if strcmp('rhd',filetype) | strcmp('mat',filetype)
          mat_fname = strcat(current_fname(1:end-4),'.mat');
          
            if ~exist(mat_fname,'file') %convert .rhd to .mat if necessary (this will slow things down the first time, but be much faster every time you want to access the .rhd in the future
                save_intan_as_mat(current_fname)
            end
            
          load(mat_fname, 'board_adc_data','frequency_parameters','t_board_adc');% faster than loading the .rhd
          songdata = board_adc_data; %entire audio trace for file
          fs = frequency_parameters.amplifier_sample_rate;
          t1_intan = t_board_adc(1); % in seconds, t=0 is when the record button was pressed

       elseif strcmp(filetype,'cbin')
          [songdata, fs] = evsoundin(cd,current_fname,'obs0'); %entire audio trace for file. assumes obs0 is audio.
       else
           error(['Something went wrong while trying to load ',current_fname,' in loop iteration ',num2str(i)])
       end
       
       % load SAUCY results - designed to be compatible with various SAUCY
       % iterations, as long as it has a properly formatted list of .txt files
       if exist('newdataobj','var')
           allspikes = transpose(newdataobj.spike_times{end}); % currently assumes last cluster contains spike.
           clear newdataobj
       elseif exist('dataobj','var')
          allspikes = dataobj.spike_times{end}; % currently assumes last cluster contains spike.
          clear dataobj
       elseif exist('Data','var')
           allspikes = Data.spiketimes{end}; % currently assumes last cluster contains spike.
           clear Data
       else
           error('Can not find SAUCY output struct')
       end
       
       time_vector = (0:numel(songdata)-1)/fs; % time vec for whole file in seconds
                   
       % find indices in time_vector corrosponding to artifact periods. 
       % create a binary vector where 1 = artifact, 0 = no artifact. if
       % mark_artificat_period was not run this is skipped
        if exist([current_fname,'.artifact_epochs.mat'])
            load([current_fname,'.artifact_epochs.mat'])
           
            % create a binary vector the same length as the file
            binary_artifact_vector = zeros(1,length(time_vector),1);
       
            for z = 1:length(artifact_offsets_in_seconds)
            artifact_idx = find(time_vector >= artifact_onsets_in_seconds(z) & time_vector <= artifact_offsets_in_seconds(z));
            binary_artifact_vector(artifact_idx) = 1;
            end
       end

       % find motif iterations within file
       seq_start = strfind(labels,syl_or_seq);
       seq_end = seq_start+(length(syl_or_seq)-1);
       pre_seq_onsets = ((onsets(seq_start))/1000) - pre_onset_time_s; % convert to seconds, add pre-syl time
       seq_onsets = ((onsets(seq_start))/1000);
       seq_offsets = ((offsets(seq_end))/1000);
       
       if isrow(pre_seq_onsets) %if labels were generated by TweetyNet/vak they are in rows, and we need them to be columns.
       pre_seq_onsets = transpose(pre_seq_onsets);
       end

       if isrow(seq_offsets)
       seq_offsets = transpose(seq_offsets);
       end
 
       %% get spikes during motif
       
       for c = 1:length(seq_onsets)
        ct = ct + 1;

        % get onsets and offsets for all syllables in motif
        all_motif_syl_onsets_idx = seq_start(c):seq_start(c)+length(syl_or_seq)-1;
        all_motif_syl_onsets = ((onsets(all_motif_syl_onsets_idx))/1000);
        
        all_motif_syl_offsets_idx = all_motif_syl_onsets_idx; % the indices are equivalent, just stored in a different array
        all_motif_syl_offsets = ((offsets(all_motif_syl_offsets_idx))/1000);
       
        %get sequence intersyllable gaps
        
        syl_gap_start = all_motif_syl_offsets(1:end-1);
        syl_gap_end = all_motif_syl_onsets(2:end);
        syl_gap_durations = syl_gap_end - syl_gap_start;
        durations_over_limit = find(syl_gap_durations > max_intersyl_s);

        
        % spikes
        spike_indices = find(allspikes >= (pre_seq_onsets(c,1)) & allspikes <= (seq_offsets(c,1) + postmotor_cushion)); % indices for spikes between seq onset and offset plus postmotor cushion
        
        % only proceed with this trial if no intersyl durations are over
        % the limit 
        
        if isempty(durations_over_limit)
            
            
       % assemble array for this syl/seq trial
        seq_spikes = allspikes(spike_indices);
        
        motif_spikes{1} = files(i,1); % originating file
        motif_spikes{2} = syl_or_seq;
        motif_spikes{3} = [all_motif_syl_onsets,all_motif_syl_offsets]; % sequence onsets and offsets for this trial
        
        recfilename(ct,1) = files(i,1);
        syl_ons_offs(ct,:) = [all_motif_syl_onsets(1),all_motif_syl_offsets(1)]; 

        neuralcase.recfilename(ct,1) = files(i,1);
        neuralcase.syl_ons_offs(ct,:) = [all_motif_syl_onsets(1),all_motif_syl_offsets(1)]; 
        neuralcase.syl_ons_offs_intan(ct,:) = [all_motif_syl_onsets(1) + t1_intan,all_motif_syl_offsets(1) + t1_intan]; 

        %assemble previous and next syllables into a string array. If there
        %is an intersection with the start or end of a file the array will
        %be empty at those positions

        prev_syls = strings(1,get_N_back_syls); %empty string array, we will populate later with syllables
        next_syls = strings(1,get_N_forward_syls); %empty string array, we will populate later with syllables
        
        if seq_start(c) - get_N_back_syls > 0 %determine if we have enough previous syllables to reach get_N_back_syls
            max_w = get_N_back_syls; % max_w is the number of previous syllables we will record
        else
            max_w = seq_start(c)-1;
        end
        
        if seq_end(c) + get_N_forward_syls < length(labels) %determine if we have enough next syllables to reach get_N_forward_syls
            max_u = get_N_forward_syls; % max_u is the number of next syllables we will record
        else
            max_u = length(labels)-seq_end(c); % seq_end(c) is the idx where our sequence ends in 'labels'
        end
        
        
        for w = 1:max_w
            prev_syls(end-w+1) =  labels(seq_start(c)-w);
        end
        
        if max_w < get_N_back_syls
            prev_syls(get_N_back_syls-max_w) = 'start file'; %indicate that data was truncated by file cutoff
            prev_syls(1:get_N_back_syls-max_w-1) = 'start file'; %get rid of empty spaces
        end
        
        for u = 1:max_u
            next_syls(u) =  labels(seq_end(c)+u);
        end
        
        if max_u < get_N_back_syls
            next_syls(max_u + 1) = 'end file'; %indicate that data was truncated by file cutoff
            next_syls(max_u + 2:end) = 'end file';
        end
        
        motif_spikes{4} = prev_syls; % next syl
        motif_spikes{5} = next_syls; % next syl
        motif_spikes{6} = seq_spikes; % entire pre-onset, seq, and post spikes with align point set to 0
               
        syls_prev(ct,:) = prev_syls;
        syls_next(ct,:) = next_syls;
        syl_in_seq(ct,:) = strcat(prev_syls(6), syl_or_seq, next_syls(1));
        spiketrains{ct} = seq_spikes;

        neuralcase.syls_prev(ct,:) = prev_syls;
        neuralcase.syls_next(ct,:) = next_syls;
        neuralcase.syl_in_seq(ct,:) = strcat(prev_syls(6), syl_or_seq, next_syls(1));
        neuralcase.spiketrains{ct,1} = seq_spikes;
        neuralcase.spiketrains_intan{ct,1} = seq_spikes + t1_intan;

        if exist('weighted_avg_labelvec','var')
            motif_spikes{7} = weighted_avg_labelvec(seq_start(c):seq_end(c));
            motif_spikes{8} = wiener_entropy_labelvec(seq_start(c):seq_end(c));
            motif_spikes{9} = amp_at_pitchquant_labelvec(seq_start(c):seq_end(c));
            motif_spikes{10} = t_assay_labelvec(seq_start(c):seq_end(c));
            motif_spikes{11} = spect_entropy_labelvec(seq_start(c):seq_end(c));

            acoustics.weighted_avg_pitch(ct,1) = weighted_avg_labelvec(seq_start(c):seq_end(c));
            acoustics.wiener_entropy(ct,1) = wiener_entropy_labelvec(seq_start(c):seq_end(c));
            acoustics.amplitude(ct,1) = amp_at_pitchquant_labelvec(seq_start(c):seq_end(c));
            acoustics.t_assay(ct,1) = t_assay_labelvec(seq_start(c):seq_end(c));
            acoustics.spect_entropy(ct,1) = spect_entropy_labelvec(seq_start(c):seq_end(c));

            neuralcase.weighted_avg_pitch(ct,1) = weighted_avg_labelvec(seq_start(c):seq_end(c));
            neuralcase.wiener_entropy(ct,1) = wiener_entropy_labelvec(seq_start(c):seq_end(c));
            neuralcase.amplitude(ct,1) = amp_at_pitchquant_labelvec(seq_start(c):seq_end(c));
            neuralcase.t_assay(ct,1) = t_assay_labelvec(seq_start(c):seq_end(c));
            neuralcase.spect_entropy(ct,1) = spect_entropy_labelvec(seq_start(c):seq_end(c));
           
        else
            motif_spikes{7} = [];
            motif_spikes{8} = [];
            motif_spikes{9} = [];
            motif_spikes{10} = [];
            motif_spikes{11} = [];
        end
        
        if exist('pitchContours','var')
            col = 12;
            motif_spikes{col} = pitchContours((seq_start(c):seq_end(c)));
            acoustics.pitch_cont(ct,1) = pitchContours((seq_start(c):seq_end(c)));
            neuralcase.pitch_cont(ct,1) = pitchContours((seq_start(c):seq_end(c)));
        else
            col = 12;
            motif_spikes{col} = [];
        end
        

       % sort into array for artifact-free trials or artifact-containing
       % trials
     if exist([current_fname,'.artifact_epochs.mat'])
        warning('still need to test that the bit code for comparing mark artifact output works correctly')
       time_minus_offset = time_vector - seq_offsets(c,1); % find closest sample to offset time
       offset_sample = find(time_minus_offset > (min(abs(time_minus_offset))-5e-06) & time_minus_offset < (min(abs(time_minus_offset)+5e-06)));
       
       time_minus_pre_time = time_vector - pre_seq_onsets(c,1);
       pre_time_onset_sample = find(time_minus_pre_time > (min(abs(time_minus_pre_time))-5e-06) & time_minus_pre_time < (min(abs(time_minus_pre_time)+5e-06)));
       
       if sum(binary_artifact_vector(pre_time_onset_sample:offset_sample)) > max_artifact_samples
         artifact_trials = vertcat(artifact_trials,motif_spikes);
       else
         clean_trials = vertcat(clean_trials,motif_spikes);
       end
    else
        clean_trials = vertcat(clean_trials,motif_spikes); % if no mark_artifact files, everything goes into clean_trials
%         clean_trials=motif_spikes;
     end
        end
       end
        

   end

%% tabulate clean_trials
if ~isempty(clean_trials)
    clean_trials_table = cell2table(clean_trials,'VariableNames',{'recname',...
        'syl_seq','seq_ons_offs','prev_syls','next_syls','seq_spikes', ...
        'pitch_wgt_avg', 'wiener_entropy','amp_at_pitchquant',...
        't_quant','spect_entropy','pitch_contours'});
else 
    clean_trials_table = table([]);
end
   
   
  %% Put into data structure and save output
if ~isempty(clean_trials_table)
    date = datetime;
    date.Format = 'yyyy-MM-dd';
    date = char(date);
    neuralcase.date_compiled = date;
   
    st = dbstack;
    generated_using_function = st.name;
    
    if ~isfolder('compile_cases')
        mkdir('compile_cases') 
    end
    
    cd('compile_cases')

       
    savename =[case_ID,'_ch',num2str(neural_channel),'_premotor_',num2str(pre_onset_time_s*1000),'ms_spiketimes_acoustics','_',date]; 
    save(savename,...
        'syl_or_seq',...
        'neural_channel',...
        'fs',...
        'pre_onset_time_s',...
        'clean_trials',...
        'artifact_trials',...
        'files',...
        'max_intersyl_s',...
        'generated_using_function',...
        'postmotor_cushion',...
        'clean_trials_table',...
        'neuralcase',...
        'date');


    cd ../
       disp(['Done. There are ',mat2str(size(clean_trials,1)),' clean trials and ',mat2str(size(artifact_trials,1)),' trials with artifact.'])
       
     toc % end timer
end

end