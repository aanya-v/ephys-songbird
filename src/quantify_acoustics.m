% This code was modified from 'headphones_quantify_pitch_lukas' 4/2013
% update. Unlike previous versions, this can take song data from either
% cbin or rhd files. It IS NOT set up for quantifying headphones
% experiments (e.g. with a minimic channel)
% 
% See 'README for headphones_quantify_pitch_lukas.doc' in
% R:\SoberLab\MATLAB_code_depository\Song_analysis\pitch_shift_analysis for
% more information on the core functionality of the code
%
%
% This function does two things:
% 
% 1.  Adds pitch (parabolic and weighted average), amplitude, spectral entropy,
%     and wiener entropy measurements to the .not.mat file accompanying each 
%    .cbin file in the current directory.
%       --> and a few other variables are saved  - see the comment block
%           "variables are now saved to .not.mat file for this song"
%
% 2.  Saves a summary file describing the pitch (parabolic and weighted average),
%     amplitude, spectral entropy, and wiener entropy of all instances of a 
%     given syllable in the directory.
%       --> and a few other variables are saved too - see the comment block
%           "organize and save data collected from all loaded .not.mat files 
%            into a single summary file"
%
% Prior to running this, you must do the following:
% 
% 1.  Collect .cbin/.rhd song files
%     
% 2.  Label syllables using evsonganaly.m or evsonganaly_intan_new.m
% 
% 3.  Write parameters for measuring the pitch of each syllable, these should be saved in a file called
%     syllable_params_by_bird_(YOUR NAME).  Use one_syl_spectrogram.m to
%     figure these parameters out. Note that you may need to add intan
%     readability, or copy Lyndie's version.
%
% 4. Set the "Variables you must set before running this function"
%    especially spbb_fname.
%
% Arguments:
%
% SYL_TO_QUANT: single syllable ('a') or syllable sequence ('abcd').  If
%   it is a single syllable, pitch/amplitude/entropy will be quantified for all
%   syllables of that type.  If it is a sequence, argument N_SYL_IN_SEQUENCE
%   determines which syllable in that particular sequence will be
%   pitch/amplitude/entropy-quantified.  Syllables not in the sequence will
%   NOT have anything quantified.
% N_SYL_IN_SEQUENCE: optional argument (default is 1).  Only set this if
%   SYL_TO_QUANT is more than 1 syllable long.  For example, if you want to
%   quantify only the second 'b' in the particular sequence 'abbcd', then
%   N_SYL_IN_SEQUENCE should be 3.  The summary file name is helpfully
%   descriptive: syl_b_number_3_in_sequence_abbcd
% FILETYPE: 'cbin' or 'rhd'
%
% Examples:
%   headphones_quantify_pitch('a',[],'rhd');
%   headphones_quantify_pitch('abbcd',3,'cbin')

function quantify_acoustics(syl_to_quant,n_syl_in_sequence,filetype)
if nargin==1, n_syl_in_sequence = 1;  end  % assume syl_to_quant is 1 syllable long or is first syllable in motif

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variables you must set before running this function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spbb_fname='syllable_params_by_bird_leila'; % Specify the name of your syllable_params_by_bird*.m file

% Important! If you are running cbin files, this code will assume that channels are set up as follows:
channel_with_song='obs0';       % cage mic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tasks before main loop: make list of songfiles & define a few variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create batchfile with all .not.mat files
if isunix
    !ls *.not.mat > batchfile
else
    !dir /B **.not.mat > batchfile
end
fid=fopen('batchfile','r');fn=fgetl(fid);fclose(fid);
if ~ischar(fn), error(['No .not.mat files found in current directory: ' pwd]); end
fid=fopen('batchfile','r'); %re-open for reading later on

% Extract date from name of first file.  This will assume that the end of
% each .cbin file is *MMDDYY_time.filenumber.cbin.
[~,day,month,year]=fn2date_rhd(fn);
date_str=[num2str(month,'%02i') '_' num2str(day,'%02i') '_' num2str(year)];

% Discovering "birdname" - assumed to be the beginning of the filename
% string up to the first underscore - e.g. all files should be named
% "birdname_xxx.cbin"
birdname=fn(1:strfind(fn,'_')-1);

% Get stuff from syllable_params_by_bird*.m --> spbb_fn
try %check if syllable_params_by_bird*.m has the fourth return value HARMONIC_NUMBER - which is used to calculate spectral entropy
    %and wiener entropy using the frequency range [PEAK_PINTERP - (half the width of fundamental frequency) : PEAK_PINTERP + (half the width of fundamental frequency)]
    %for each syllable - where fundamental frequency = PEAK_PINTERP/HARMONIC_NUMBER.
    eval(sprintf('[f_cutoff,t_assay_original,spect_params,harmonic_number]=%s(birdname,syl_to_quant(n_syl_in_sequence));',spbb_fname))
catch e %if not, just load the three other return values, and use F_CUTOFF as the frequency range to calculate spectral entropy and wiener entropy.
    warning(['Function ' spbb_fname ' did not return HARMONIC_NUMBER for bird ' birdname ', syllable ' syl_to_quant(n_syl_in_sequence),...
        '. SPECT_ENTROPY and WIENER_ENTROPY will be calculated using frequency range [F_CUTOFF(1) : F_CUTOFF(2)] ',...
        'instead of [PEAK_PINTERP - (half the width of fundamental frequency) : PEAK_PINTERP + (half the width of fundamental frequency)].']);
    eval(sprintf('[f_cutoff,t_assay_original,spect_params]=%s(birdname,syl_to_quant(n_syl_in_sequence));',spbb_fname))
end
if strcmp(f_cutoff,'undefined') || strcmp(t_assay_original,'undefined') || strcmp(spect_params,'undefined')
    error(['Syllable ' syl_to_quant(n_syl_in_sequence) ' for bird ' birdname ' does not have all ',...
        'three required parameters F_CUTOFF, T_ASSAY, SPECT_PARAMS defined in function ' spbb_fname]);
end

%Talkback
disp(['This is for bird ' birdname])
disp(['Window duration = ' num2str(spect_params(2)) ' msec, overlap = ' num2str(spect_params(1)),...
    ', assaying at t=' num2str(t_assay_original) ' with freq cutoff at [' num2str(f_cutoff) ']'])

        f_lo = 500; f_hi = 10000;
        f_cutoff_full = [f_lo f_hi];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main loop: load one file at a time, save some outputs back to the
% original .not.mat file, and keep others to be saved in a single file
% containing info from all .not.mat files in the current directory.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ct=1; %count for saving to summary file; incremented each time a syllable is processed
t_out = [];
while (1)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Before calculating syllable-by-syllable stats, load the file's raw 
    % data and existing .not.mat data.  Check validity of a couple things.
    % Obtain information about target syllable ids, syllable onsets, 
    % offsets, song times. Obtain pitch quantification parameters.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Current file names.
    fn=fgetl(fid);                          % get name of a .not.mat filename
    cbin_fn=fn(1:end-8);                    % name of corresponding .cbin/rhd file
    rec_fn=[cbin_fn(1:end-5) '.rec'];       % name of corresponding .rec file
    if (~ischar(fn));fclose(fid);break;end  % Break while loop when last filename is reached
    
    % if file was recorded on a different day than the first file in directory,
    % summary file's date is 'multiple days'.
    [~,day,month,year]=fn2date_rhd(fn);
    date_str_current_file=[num2str(month,'%02i') '_' num2str(day,'%02i') '_' num2str(year)];
    if ~strcmp(date_str_current_file,date_str)
        date_str = 'multiple_days';
    end
    
    % Load .not.mat file with specified fields so they are resaved without overwriting data from other syllables
    % Onsets, offsets, labels are not resaved later because nothing will ever change for them in this function.
    warning off % suppress warnings if vars are not found in file because
                % this is expected the first time headphones_quantify_pitch
                % is run on a .not.mat file.
    eval(sprintf(['load(''%s'',''onsets'',''offsets'',''labels'',',...
    '''P1_save_allfreqs_labelvec'',''F1_save_allfreqs_labelvec'',''P1_save_labelvec'',''F1_save_labelvec'',',...
    '''peak_pinterp_labelvec'',''peak_pinterp_labelvec_shifted'',''peak_pinterp_labelvec_minimic'',',...
    '''weighted_avg_labelvec'',''weighted_avg_labelvec_shifted'',''weighted_avg_labelvec_minimic'',',...
    '''amp_at_pitchquant_labelvec'',''amp_at_pitchquant_labelvec_shifted'',''amp_at_pitchquant_labelvec_minimic'',',...
    '''spect_entropy_labelvec'',''spect_entropy_labelvec_shifted'',''spect_entropy_labelvec_minimic'',',...
    '''wiener_entropy_labelvec'',''wiener_entropy_labelvec_shifted'',''wiener_entropy_labelvec_minimic'',',...
    '''t_assay_labelvec'',''pitchContours'')'],fn));
    warning on
    
    %Added 6/28/2013. This fixes a bug which happens in the following sequence.
    %(1) Run headphones_quantify_pitch_lukas once
    %(2) You decide to relabel some syllables in evsonganaly or change threshold
    %(3) This change alters the number of labels in some files.
    %       if you relabeled some syllables this would occur if you used the
    %       'Edit' button in evonsongaly to split or combine syllables.
    %       if you changed syllable thresholds the function could
    %       automatically split or combine syllables with the new threshold.
    %(4) Run headphones_quantify_pitch_lukas again.
    %(5) At this spot in the code, the loaded variables (such as
    %weighted_avg_labelvec) DO NOT HAVE THE SAME # OF ELEMENTS AS THE
    %LABELS STRING!!  THIS CAUSES THE PITCHES TO BE PLACED IN THE WRONG
    %INDEX AND SAVED BACK TO FILE!!
    %
    %The code below fixes this bug by ensuring that every loaded variable
    %that was previously saved by headphones_quantify_pitch_lukas has
    %the same number of elements as the labels vector.  If not, it clears
    %that variable and the rest of the code proceeds as if
    %headphones_quantify_pitch is being run for the first time.  Then the
    %correct variables will be saved.
    lblsize = size(labels,2);
    if exist('P1_save_allfreqs_labelvec','var') && size(P1_save_allfreqs_labelvec,2)~=lblsize, clear P1_save_allfreqs_labelvec; end
    if exist('F1_save_allfreqs_labelvec','var') && size(F1_save_allfreqs_labelvec,2)~=lblsize, clear F1_save_allfreqs_labelvec; end
    if exist('P1_save_labelvec','var') && size(P1_save_labelvec,2)~=lblsize, clear P1_save_labelvec; end
    if exist('F1_save_labelvec','var') && size(F1_save_labelvec,2)~=lblsize, clear F1_save_labelvec; end
    if exist('peak_pinterp_labelvec','var') && size(peak_pinterp_labelvec,2)~=lblsize, clear peak_pinterp_labelvec; end
    if exist('peak_pinterp_labelvec_shifted','var') && size(peak_pinterp_labelvec_shifted,2)~=lblsize, clear peak_pinterp_labelvec_shifted; end
    if exist('peak_pinterp_labelvec_minimic','var') && size(peak_pinterp_labelvec_minimic,2)~=lblsize, clear peak_pinterp_labelvec_minimic; end    
    if exist('weighted_avg_labelvec','var') && size(weighted_avg_labelvec,2)~=lblsize, clear weighted_avg_labelvec; end
    if exist('weighted_avg_labelvec_shifted','var') && size(weighted_avg_labelvec_shifted,2)~=lblsize, clear weighted_avg_labelvec_shifted; end
    if exist('weighted_avg_labelvec_minimic','var') && size(weighted_avg_labelvec_minimic,2)~=lblsize, clear weighted_avg_labelvec_minimic; end
    if exist('amp_at_pitchquant_labelvec','var') && size(amp_at_pitchquant_labelvec,2)~=lblsize, clear amp_at_pitchquant_labelvec; end
    if exist('amp_at_pitchquant_labelvec_shifted','var') && size(amp_at_pitchquant_labelvec_shifted,2)~=lblsize, clear amp_at_pitchquant_labelvec_shifted; end
    if exist('amp_at_pitchquant_labelvec_minimic','var') && size(amp_at_pitchquant_labelvec_minimic,2)~=lblsize, clear amp_at_pitchquant_labelvec_minimic; end    
    if exist('spect_entropy_labelvec','var') && size(spect_entropy_labelvec,2)~=lblsize, clear spect_entropy_labelvec; end
    if exist('spect_entropy_labelvec_shifted','var') && size(spect_entropy_labelvec_shifted,2)~=lblsize, clear spect_entropy_labelvec_shifted; end    
    if exist('spect_entropy_labelvec_minimic','var') && size(spect_entropy_labelvec_minimic,2)~=lblsize, clear spect_entropy_labelvec_minimic; end   
    if exist('wiener_entropy_labelvec','var') && size(wiener_entropy_labelvec,2)~=lblsize, clear wiener_entropy_labelvec; end      
    if exist('wiener_entropy_labelvec_shifted','var') && size(wiener_entropy_labelvec_shifted,2)~=lblsize, clear wiener_entropy_labelvec_shifted; end       
    if exist('wiener_entropy_labelvec_minimic','var') && size(wiener_entropy_labelvec_minimic,2)~=lblsize, clear wiener_entropy_labelvec_minimic; end  
    if exist('t_assay_labelvec','var') && size(t_assay_labelvec,2)~=lblsize, clear t_assay_labelvec; end          

    % Load raw song waveforms
    if contains(filetype,'cbin')
    [rawsong, Fs]=evsoundin('.', cbin_fn,channel_with_song);
    elseif contains(filetype,'rhd')
        matname = [cbin_fn(1:end-3),'mat'];
        load(matname)
        rawsong = board_adc_data;
        Fs = frequency_parameters.amplifier_sample_rate;
    elseif contains(filetype,'mat')
        matname = [cbin_fn(1:end-3),'mat'];
        load(matname)
        rawsong = board_adc_data;
        Fs = frequency_parameters.amplifier_sample_rate;
    else
        error('invalid filetype')
    end
    
    %Obtain syllable onsets, offsets, and time of each value in rawsong   
    onsets=onsets/1000;offsets=offsets/1000;% get onsets and offsets into sec, not msec
    syl_dur = offsets-onsets;
    t=[1:length(rawsong)]/Fs; % Time
    
    % Find (id) - index of target syllables.
    % If this targeted syllable has additional restrictions (e.g. look for
    % a syllable only within a certain sequence), only syllables within 
    % this sequence are quantified.  Summary file gets an extra string at
    % the end to indicate which motif and syllable.
    id = strfind(labels,syl_to_quant);  %find either single syllable 'a' or motif 'aabbc'
    id = id + n_syl_in_sequence - 1; %change ID to match # syllable in sequence 'aabbcc', if single syllable 'a' ID stays the same.
    disp(['Found ' num2str(length(id)) ' examples of ''' syl_to_quant ''' in ' fn])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Individual-syllable loop: for each targeted syllable in the file
    % indexed by id: calculate pitch (parabolic interpolation), pitch
    % (weighted average), amplitude, spectral entropy, wiener entropy,
    % frequency slice used in some calculations (F1_save), power at those
    % frequency slices (P1_save), non-cutoff frequencies and powers
    % (F1/P1_save_allfreqs), shift amount taken from filename, precise
    % time of recording (t_out).
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ct_within_each_file=1;  % count for saving to each .not.mat file, as opposed to total count (ct) for summary file

    onset_within_file = zeros(length(id));
    offset_within_file = zeros(length(id));
    syl_duration_within_file = zeros(length(id),1);
    syl_wav_within_file = cell(length(id),1);    

    for x=1:length(id) %for each targeted syllable in file

        % Define onsets, offsets, and time of pitch/amplitude quantification
        on=onsets(id(x)); off=offsets(id(x));
        if t_assay_original>1   % If t_assay>1, then it is used as a % of the duration of the syl
            t_assay=.01*t_assay_original*(off-on);
        else % Otherwise it is used as seconds since start of syllable
            t_assay=t_assay_original;
        end

        t_assay_within_file(ct_within_each_file) = t_assay; %this vector saves time in seconds since start of syllable for all quantified syllables
        
        if (off-on)<spect_params(2)/1000    % if syllable is too short for pitch analysis
            off=on+spect_params(2)/1000;    % extend end of syllable
            disp('Extending offset')
        end
        
        onset(ct) = on;
        onset_within_file(ct_within_each_file) = on;
        offset(ct) = off;
        offset_within_file(ct_within_each_file) = off;
        syl_duration(ct) = off-on;
        syl_duration_within_file(ct_within_each_file) = off-on;
        
        %Obtain raw waveform of syllable using onset/offset times
        on_id=find(abs((t-on))==min(abs(t-on))); % what time in rawsong syllable started
        off_id=find(abs((t-off))==min(abs(t-off))); % ... and ended
        syl_wav=rawsong(on_id:off_id); % waveform of the syllable
        
        syllable_wav{ct} = syl_wav;
        syl_wav_within_file{ct_within_each_file} = syl_wav;
        
        % Compute spectrogram for whole syllable
        [S1,F1,T1,P1] =spect_from_waveform(syl_wav,Fs,0,spect_params);
        if ~exist('save_T1','var'), save_T1=[]; end
        if length(T1)>length(save_T1) % save_T1 is the longest T1
            save_T1=T1; %times at which spectrogram was computed
        end
        
        % In the array of times that fft (fast fourier transform) was
        % taken(T1), find which time is closest to time to assay (t_assay).
        t_id=find(abs((T1-t_assay))==min(abs(T1-t_assay)));
%         t_assay_se = [t_assay-0.016 t_assay
%         t_id_sp_en = t_id-16:1:t_id+16;
        if length(t_id)>1, t_id=min(t_id); end
        
        % Save the slice of the whole spectrogram at time of pitch quantification
        P1_save_allfreqs(ct,:)=P1(:,t_id)'; % for summary file, each row is spectrogram slice's power for 1 syllable
        P1_save_allfreqs_within_file(ct_within_each_file,:)=P1(:,t_id)'; % for .not.mat file for this song
        F1_save_allfreqs(ct,:)=F1'; % for summary file, each row is spectrogram slice's frequencies for above
        F1_save_allfreqs_within_file(ct_within_each_file,:)=F1'; % for .not.mat file for this song
        
        % Reduce F1 and P1 to the range of frequencies specified for that
        % syllable by syllable_parameters_by_bird_xx.m
        f_cut_id=find(F1>f_cutoff(1) & F1<f_cutoff(2));
        F1_for_pitchquant=F1(f_cut_id);
        P1_for_pitchquant=P1(f_cut_id,:);
        spect_slice=P1_for_pitchquant(:,t_id);
        
        % Quantify SECOND harmonic
%         f_cut_id_2=find(F1>f_cutoff(3) & F1<f_cutoff(4));
%         F1_for_pitchquant_2=F1(f_cut_id_2);
%         P1_for_pitchquant_2=P1(f_cut_id_2,:);
%         spect_slice_2=P1_for_pitchquant_2(:,t_id);
        
        % spect for the fuller frequency range of the syllable
        f_lo = 500; f_hi = 10000;
        f_cutoff_full = [f_lo f_hi];
        f_all_id = find(F1>f_lo & F1<f_hi);
        F1_for_quant=F1(f_all_id);
        P1_for_quant=P1(f_all_id,:);
        spect_all = P1_for_quant(:,t_id)';
        
        % Save the region of spectrum used to compute pitch at time of pitch quantification
        P1_save(ct,:)=spect_slice; % for summary file, each row is spectrogram at time of pitch quantification over the frequencies specified by syllable_parameters_by_bird_xx.m
        P1_save_within_file(ct_within_each_file,:)=spect_slice;
        F1_save(ct,:)=F1_for_pitchquant; % for summary file, each row is spectrogram slice's frequencies for above
        F1_save_within_file(ct_within_each_file,:)=F1_for_pitchquant;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Computation of pitch/amplitude/entropy of individual
        % syllable, for the bird's actual song (channel_with_song).
        % For all calculations except amplitude, we use the spect_slice
        % obtained above - constrained frequency range in a single time bin.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        col = 1; %store in column 1 because channel # 1
        
        % Internal function compute_scalar_pitch extracts pitch using
        % parabolic interpolation around convex peak
        pitch_tmp=compute_scalar_pitch(spect_slice,F1_for_pitchquant);
        peak_pinterp(ct,col) = pitch_tmp; % for summary file, column 1 is pitches for syllables in channel_with_song
        peak_pinterp_within_file(ct_within_each_file,col) = pitch_tmp; %for .not.mat file for this song
        
        % Second Harmonic
%         pitch_tmp_2=compute_scalar_pitch(spect_slice_2,F1_for_pitchquant_2);
%         peak_pinterp_2(ct,col) = pitch_tmp_2; % for summary file, column 1 is pitches for syllables in channel_with_song
%         peak_pinterp_within_file_2(ct_within_each_file,col) = pitch_tmp_2; %for .not.mat file for this song
        
        % Internal function compute_scalar_pitch_weighted_average extracts
        % pitch using weighted average of frequencies around convex peak.
        pitch_tmp_weighted_avg=compute_scalar_pitch_weighted_avg(spect_slice,F1_for_pitchquant);
        weighted_avg(ct,col) = pitch_tmp_weighted_avg; % for summary file, column 1 is pitches for syllables in channel_with_song
        weighted_avg_within_file(ct_within_each_file,col) = pitch_tmp_weighted_avg; %for .not.mat file for this song
        
        % Second Harmonic
%         pitch_tmp_weighted_avg_2=compute_scalar_pitch_weighted_avg(spect_slice_2,F1_for_pitchquant_2);
%         weighted_avg_2(ct,col) = pitch_tmp_weighted_avg_2; % for summary file, column 1 is pitches for syllables in channel_with_song
%         weighted_avg_within_file_2(ct_within_each_file,col) = pitch_tmp_weighted_avg_2; %for .not.mat file for this song

        % Smooth sound waveform and compute amplitude with internal
        % function compute_summed_amp_window
        sm_win = 15.0; % ms
        sm=evsmooth(syl_wav,Fs,0.01,512,0.8,sm_win,f_lo,f_hi);
        [amp_tmp, t_assay_amp]=compute_summed_amp_window(sm,t_assay,Fs);  % Compute summed amplitude in 16msec around pitch quant
        amp_at_pitchquant(ct,col)=amp_tmp; % for summary file: this will continue to be added to as multiple files are loaded
        amp_at_pitchquant_within_file(ct_within_each_file,col)=amp_tmp; % for .not.mat file for current song: this is cleared each time a new file is loaded
%         if amp_tmp > 0.25
%             whos sm amp_tmp;figure;plot(sm); hold on; plot(amp_tmp); return   
%         end
        % Obtain spectral slice for calculating spectral and Wiener entropy
        if exist('harmonic_number','var')
            % For calculating spectral and Wiener entropy, obtain a SPECT_SLICE with frequency range
            % [PEAK_PINTERP - (half the width of fundamental frequency) : PEAK_PINTERP + (half the width of fundamental frequency)],
            % where fundamental frequency = PEAK_PINTERP/HARMONIC_NUMBER. This is the way it was done in Sam Sober's 2010 paper
            f_1 = pitch_tmp/harmonic_number; %f1 = fundamental frequency for this syllable
            range = pitch_tmp + [-f_1/2 +f_1/2]; %a frequency window centered on PEAK_PINTERP, going halfway down to the next-lowest harmonic and halfway up to the next-highest harmonic
            f_cut_id_entropy = find(F1>range(1) & F1<range(2));
            P1_for_pitchquant_entropy = P1(f_cut_id_entropy,:);
            spect_slice_entropy = P1_for_pitchquant_entropy(:,t_id);
            abs_S1 = abs(spect_slice_entropy);
%         else %if syllable_params_by_bird didn't return HARMONIC_NUMBER, use SPECT_SLICE, which has F_CUTOFF as the frequency range       
%             abs_S1 = abs(spect_slice); 
        else %if syllable_params_by_bird didn't return HARMONIC_NUMBER, use SPECT_ALL, which has frequency range between 0 and 8kHz (2010 Wolgemuth Sober paper)       
            abs_S1 = abs(spect_all);
        end
        
        % Compute spectral entropy
        abs_S1_rat = abs_S1/sum(abs_S1);
        spect_entropy_tmp = -sum(abs_S1_rat.*log10(abs_S1_rat));
        spect_entropy(ct,col)=spect_entropy_tmp; % for summary file, column 1 is entropy for syllables in channel_with_song
        spect_entropy_within_file(ct_within_each_file,col)=spect_entropy_tmp; %for .not.mat file for this song
        
        % Compute Wiener entropy
        wiener_entropy_tmp=10*log10(geomean(abs_S1)/mean(abs_S1)); %in dB
        wiener_entropy(ct,col)=wiener_entropy_tmp; % for summary file, column 1 is entropy for syllables in channel_with_song
        wiener_entropy_within_file(ct_within_each_file,col)=wiener_entropy_tmp; %for .not.mat file for this song

        % Extract Pitch Contours - derived from external function
        % pitch_contours{ct} = pitchContours(id(x));
        % pitch_contours_within_file{ct_within_each_file} = pitchContours(id(x));
       
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Additional things to save to summary file: filename for current
        % syllable, recording time (t_out), shift amount.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             
        %Filename for current syllable
        fname_arr{ct}=cbin_fn;

        % Compute time of recording - FROM .REC HEADER
        if contains(cbin_fn,'rhd')
            yearmonthday = strtok(cbin_fn(end-16:end),'_');
            hrminsec = strtok(cbin_fn(end-9:end),'.');
            timestamp = [yearmonthday,hrminsec,num2str(floor(1000*(on-floor(on))),'%03i')];
            t_out(ct)=str2num(timestamp); 
        elseif contains(cbin_fn,'mat')
            idx = strfind(cbin_fn,'_'); 
            yearmonthday = cbin_fn(idx(1)+1:idx(2)-1); %after first "_" is the date
            hrminsec = fn(idx(2)+1:idx(3)-1);
%             yearmonthday = strtok(cbin_fn(end-16:end),'_');
%             hrminsec = strtok(cbin_fn(end-9:end),'.');
            timestamp = [yearmonthday,hrminsec,num2str(floor(1000*(on-floor(on))),'%03i')];
            t_out(ct)=str2num(timestamp); 
        else
            
        recdata=readrecf(rec_fn);
        t_header_str=recdata.header{1};
        space_id=find(t_header_str==' ');
        t_str=t_header_str((space_id(end)+1):end);
        hr=str2num(t_str(1:2));
        minute=str2num(t_str(4:5));
        second=str2num(t_str(7:8));
        
        %the exact time down to the millisecond when the syllable occurred.
        %T_OUT format: YYYYMMDDHHmmSSSXXX Y=year, M=month, D=day, H=hour, m=minute, SSS = second (3 digits in case it was a very long file >100 seconds), XXX=milliseconds
        timestamp = [num2str(year) num2str(month,'%02i') num2str(day,'%02i'),... %also from the .rec file
            num2str(hr,'%02i') num2str(minute,'%02i') num2str(second + floor(on),'%03i'),... %3 digits for the second in case it is a very long file > 100 seconds
            num2str(floor(1000*(on-floor(on))),'%03i')]; %this line finds # of milliseconds
        timestamp = str2num(['uint64(' timestamp ')']); %convert to uint64 because it can store exact #s up to 20 digits. timestamp has 18 digits
        t_out(ct)=timestamp; %used later to sort summary file data by time
        end
        
        % TABULATE
        syltable{ct,1} = cbin_fn;
        syltable{ct,2} = id(x);
        syltable{ct,3} = syl_to_quant;
        syltable{ct,4} = [on off];
        syltable{ct,5} = syl_duration_within_file(x);
        syltable{ct,6} = f_cutoff;
        syltable{ct,7} = f_cutoff_full;
        syltable{ct,8} = sm_win;
        syltable{ct,9} = spect_params;
        syltable{ct,10} = t_assay;
%         syltable{ct,9} = save_T1;
%         syltable{ct,10} = P1_save_allfreqs_within_file(x,:);
%         syltable{ct,11} = F1_save_allfreqs_within_file(x,:);
%         syltable{ct,12} = P1_save_within_file(x,:);
%         syltable{ct,13} = F1_save_within_file(x,:);
        syltable{ct,11} = peak_pinterp_within_file(x);
        syltable{ct,12} = weighted_avg_within_file(x);
        syltable{ct,13} = amp_at_pitchquant_within_file(x);
        syltable{ct,14} = spect_entropy_within_file(x);
        syltable{ct,15} = wiener_entropy_within_file(x);
        % syltable{ct,16} = pitch_contours_within_file(x);
        
        ct=ct+1;
        ct_within_each_file=ct_within_each_file+1;
    end  % end of loop going through targeted syllables for 1 file
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Variables are now saved to .not.mat file for this song:
    % channel_with_song
    % P1_save_allfreqs_labelvec (only for channel_with_song)
    % F1_save_allfreqs_labelvec (only for channel_with_song)
    % P1_save_labelvec (only for channel_with_song)
    % F1_save_labelvec (only for channel_with_song)
    % peak_pinterp_labelvec (and _shifted and _minimic)
    % weighted_avg_labelvec (and _shifted and _minimic) 
    % amp_at_pitchquant_labelvec (and _shifted and _minimic)
    % spect_entropy_labelvec (and _shifted and _minimic)
    % wiener_entropy_labelvec (and _shifted and _minimic)
    % t_assay_labelvec
    % f_cutoff_syl_x (where x is syl_to_quant(n_syl_in_sequence))
    % t_assay_syl_x (where x is syl_to_quant(n_syl_in_sequence))
    % spect_params_syl_x (where x is syl_to_quant(n_syl_in_sequence))
    %
    % * If these previously existed in the .not.mat file, stored values will
    %   not be overwritten except those corresponding to syl_to_quant.  So if
    %   you run headphones_quantify_pitch('a') and then
    %   headphones_quantify_pitch('b'), then for example the
    %   peak_pinterp_labelvec values for 'a' will not be overwritten.  
    % * However, if you call headphones_quantify_pitch('b') twice, the values
    %   corresponding to syllable 'b' will be overwritten.
    % * If you call headphones_quantify_pitch('abb',2) and then
    %   headphones_quantify_pitch('abb',3), then the first b in motif abb
    %   will not be overwritten, because the two function calls are 
    %   quantifying syllables at two distinct positions within the motif.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Note the code below will still work if no targeted syllables found in
    % file (numel(id)==0) because it will just save the same vectors back,
    % or save new empty vectors.

    %Make new vectors if they were not loaded from .not.mat file
    if ~exist('P1_save_allfreqs_labelvec','var'), P1_save_allfreqs_labelvec = cell(size(labels)); end
    if ~exist('F1_save_allfreqs_labelvec','var'), F1_save_allfreqs_labelvec = cell(size(labels)); end
    if ~exist('P1_save_labelvec','var'), P1_save_labelvec = cell(size(labels)); end
    if ~exist('F1_save_labelvec','var'), F1_save_labelvec = cell(size(labels)); end
    if ~exist('peak_pinterp_labelvec','var'), peak_pinterp_labelvec=zeros(size(labels)); end
    if ~exist('weighted_avg_labelvec','var'), weighted_avg_labelvec=zeros(size(labels)); end
    if ~exist('amp_at_pitchquant_labelvec','var'), amp_at_pitchquant_labelvec=zeros(size(labels)); end
    if ~exist('spect_entropy_labelvec','var'), spect_entropy_labelvec=zeros(size(labels)); end
    if ~exist('wiener_entropy_labelvec','var'), wiener_entropy_labelvec=zeros(size(labels)); end
    if ~exist('t_assay_labelvec','var'), t_assay_labelvec=zeros(size(labels)); end
    if ~exist('syl_duration_labelvec','var'), syl_duration_labelvec=zeros(size(labels));end
    if ~exist('onset_labelvec','var'), onset_labelvec=zeros(size(labels)); end
    if ~exist('offset_labelvec','var'), offset_labelvec=zeros(size(labels));end

    
    % Fill in vectors with "within_file" values.  Other syllables are not
    % overwritten if they already had values in those vectors.
    for x=1:length(id)
        P1_save_allfreqs_labelvec{id(x)} = P1_save_allfreqs_within_file(x,:);
        F1_save_allfreqs_labelvec{id(x)} = F1_save_allfreqs_within_file(x,:);
        P1_save_labelvec{id(x)} = P1_save_within_file(x,:);
        F1_save_labelvec{id(x)} = F1_save_within_file(x,:);
        col = 1;
        peak_pinterp_labelvec(id(x)) = peak_pinterp_within_file(x,col);
        weighted_avg_labelvec(id(x)) = weighted_avg_within_file(x,col);
        amp_at_pitchquant_labelvec(id(x)) = amp_at_pitchquant_within_file(x,col);
        spect_entropy_labelvec(id(x)) = spect_entropy_within_file(x,col);
        wiener_entropy_labelvec(id(x)) = wiener_entropy_within_file(x,col);
        t_assay_labelvec(id(x)) = t_assay_within_file(x);
        onset_labelvec(id(x)) = onset_within_file(x);
        offset_labelvec(id(x)) = offset_within_file(x);
        syl_duration_labelvec(id(x)) = syl_duration_within_file(x);
        syl_wav_labelvec{id(x)} = syl_wav_within_file{x};
    end

    %Used to save syllable_params_by_bird parameters for this syllable to .not.mat files
    syl = syl_to_quant(n_syl_in_sequence);
    eval(sprintf('f_cutoff_pitch_syl_%s = f_cutoff;', syl)); %frequency range for pitch & pitch contours quant
    eval(sprintf('f_cutoff_syl_%s = f_cutoff_full;', syl)); %frequency range for spect/wiener entropy & amp quant
    eval(sprintf('t_assay_syl_%s = t_assay_original;', syl));
    % eval(sprintf('t_assay_amp_syl_%s = t_assay_amp;', syl));
    eval(sprintf('spect_params_syl_%s = spect_params;', syl));
    eval(sprintf('syl_dur_%s = syl_dur;',syl));

    % Save .not.mat file with the new fields, not overwriting the other
    % fields (the -append flag) or existing values within fields from other
    % syllables (because the same variables were loaded from .not.mat file
    % at top of loop).
%     eval(sprintf(['save(''%s'',''channel_with_song'',',...
%         '''P1_save_allfreqs_labelvec'',''F1_save_allfreqs_labelvec'',''P1_save_labelvec'',''F1_save_labelvec'',',...
%         '''peak_pinterp_labelvec'',''weighted_avg_labelvec'',''amp_at_pitchquant_labelvec'',',...
%         '''spect_entropy_labelvec'',''wiener_entropy_labelvec'',''t_assay_labelvec'',',...
%         '''f_cutoff_pitch_syl_%s'',''f_cutoff_syl_%s'', ''t_assay_syl_%s'', ''t_assay_amp_syl_%s'', ''spect_params_syl_%s'', ''syl_dur_%s'',''-append'')'],fn, syl, syl, syl, syl,syl,syl));
    
%     eval(sprintf(['save(''%s'',''channel_with_song'',',...
%         '''P1_save_allfreqs_labelvec'',''F1_save_allfreqs_labelvec'',''P1_save_labelvec'',''F1_save_labelvec'',',...
%         '''peak_pinterp_labelvec'',''weighted_avg_labelvec'',''amp_at_pitchquant_labelvec'',',...
%         '''spect_entropy_labelvec'',''wiener_entropy_labelvec'',''t_assay_labelvec'',',...
%         '''f_cutoff_pitch_syl_%s'',''f_cutoff_syl_%s'', ''t_assay_syl_%s'', ''spect_params_syl_%s'',',...
%         'syl_dur_%s'',''-append'')'],fn, syl, syl, syl, syl,syl));
    eval(sprintf(strcat("save('%s', 'channel_with_song','P1_save_allfreqs_labelvec',",...
        "'F1_save_allfreqs_labelvec','P1_save_labelvec','F1_save_labelvec',",...
        "'peak_pinterp_labelvec','weighted_avg_labelvec','amp_at_pitchquant_labelvec',",...
        "'spect_entropy_labelvec','wiener_entropy_labelvec','t_assay_labelvec',",...
        "'f_cutoff_pitch_syl_%s','f_cutoff_syl_%s','t_assay_syl_%s',",...
        "'spect_params_syl_%s','syl_dur_%s','-append')"),fn, syl,syl,syl, syl, syl));

    % Cleaning up so that new variables don't carry over when next file is added
    clear P1_save_allfreqs_labelvec F1_save_allfreqs_labelvec P1_save_labelvec F1_save_labelvec
    clear peak_pinterp_labelvec peak_pinterp_labelvec_shifted peak_pinterp_labelvec_minimic
    clear weighted_avg_labelvec weighted_avg_labelvec_shifted weighted_avg_labelvec_minimic
    clear amp_at_pitchquant_labelvec amp_at_pitchquant_labelvec_shifted amp_at_pitchquant_labelvec_minimic
    clear spect_entropy_labelvec spect_entropy_labelvec_shifted spect_entropy_labelvec_minimic
    clear wiener_entropy_labelvec wiener_entropy_labelvec_shifted wiener_entropy_labelvec_minimic
    clear t_assay_labelvec 
    clear P1_save_allfreqs_within_file F1_save_allfreqs_within_file P1_save_within_file F1_save_within_file
    clear peak_pinterp_within_file weighted_avg_within_file amp_at_pitchquant_within_file spect_entropy_within_file wiener_entropy_within_file
    clear t_assay_within_file 
    clear onset_labelvec onset_within_file offset_labelvec offset_within_file
    clear syl_duration_labelvec syl_duration_within_file
    clear syl_wav_labelvec syl_wav_within_file

end %end of main loop going through all files

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Organize and save data collected from all loaded .not.mat files into a
% single summary file for the specified syl_to_quant (summary_fn):
%
% P1_save_allfreqs (only for channel_with_song)
% F1_save_allfreqs (only for channel_with_song)
% P1_save (only for channel_with_song)
% F1_save (only for channel_with_song)
% peak_pinterp (column 1 = channel_with_song, 2 = _shifted, 3 =_mini_mic)
% weighted_avg (column 1 2 3 as above)
% amp_at_pitchquant (column 1 2 3 as above)
% spect_entropy (column 1 2 3 as above)
% wiener_entropy (column 1 2 3 as above)
% fname_arr
% t_out (format: YYYYMMDDHHmmSSSXXX Y=year, M=month, D=day, H=hour, m=minute, SSS = second (3 digits in case it was a very long file >100 seconds), XXX=milliseconds)
% shift
% save_T1
% t_offset_from_xcorr
% t_offset_from_xcorr_minimic
% spect_params
% f_cutoff
% t_assay
% birdname
%
% * Any existing summary file for this syllable will be overwritten!
%
% * Warning: if you quantified an uppercase syllable and re-ran this function
%   on a lowercase syllable (such as 'A', then 'a'), you will notice there is
%   only one summary file *syl_A*.  If you run this function with 'a' then
%   'A' there is only one summary file *syl_a*.  This is an unavoidable bug
%   in Windows operating system.  Recommendation: run a script to
%   automatically relabel all of your uppercase syllables to a lowercase
%   syllable that is not used to label any other syllable ('A' --> 'z').
%   "Because the Windows operating system considers two files with the same name
%   to be the same file (regardless of case), you cannot have two files with the 
%   same name in the same folder. If you save MYFILE, and myfile.mat already 
%   exists in the current folder, then MYFILE.MAT replaces myfile.mat without warning.
%   If you save myfile, and MYFILE.mat already exists in the current folder, the
%   contents of myfile.mat replace the contents of MYFILE.mat, but the name remains MYFILE.mat."
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Make name of summary file
%Note: date_str is either MM_DD_YYYY or 'multiple_days'
summary_fn = [birdname '_' date_str '_syl_' syl_to_quant(n_syl_in_sequence)];
if length(syl_to_quant)>1, summary_fn = [summary_fn '_number_' num2str(n_syl_in_sequence) '_in_sequence_' syl_to_quant]; end
disp(['Saving summary file: ' summary_fn '.mat']);

syl_sum = cell2table(syltable,'VariableNames',{'recname','label_idx',...
    'syllable','onset_offset','syl_dur','f_cutoff','f_cutoff_full',...
    'amp_sm_window','spect_params','t_assay','peak_pinterp','weighted_avg','amp_at_pitchquant',...
    'spect_entropy','wiener_entropy'});

% Sort all variables by t_out, so that for example peak_pinterp(1,:)
% corresponds to the earliest targeted syllable in directory,
% peak_pinterp(end,:) is the latest, in terms of year/month/day/hour/minute/second/millisecond sung.
if ct>1 % at least one targeted syllable was found in entire directory
    [~,t_id] = sort(t_out);

    P1_save_allfreqs = P1_save_allfreqs(t_id,:);
    F1_save_allfreqs = F1_save_allfreqs(t_id,:);
    P1_save = P1_save(t_id,:);
    F1_save = F1_save(t_id,:);
    peak_pinterp = peak_pinterp(t_id,:);
    weighted_avg = weighted_avg(t_id,:);
    amp_at_pitchquant = amp_at_pitchquant(t_id,:);    
    spect_entropy = spect_entropy(t_id,:);
    wiener_entropy = wiener_entropy(t_id,:);
    fname_arr = fname_arr(t_id);
    t_out = t_out(t_id);
    syl_duration = syl_duration(t_id);
    %save_T1 doesn't need to be sorted
    % Note that spect_params, f_cutoff, t_assay, birdname were
    % already defined earlier in function.
else %no targeted syllables found at all - make empty variables
    P1_save_allfreqs = [];
    F1_save_allfreqs = [];
    P1_save = [];
    F1_save = [];
    peak_pinterp = [];
    weighted_avg = [];
    amp_at_pitchquant = [];    
    spect_entropy = [];
    wiener_entropy = [];
    fname_arr = {};
    t_out = [];
    save_T1 = [];
    % Note that spect_params, f_cutoff, t_assay_original, birdname were
    % already defined earlier in function.
end
t_assay = t_assay_original; %save the t_assay variable as the original t_assay number that was loaded from spbb_filename

% Save summary file for the quantified syllable

% Warning: if you quantified an uppercase syllable and re-ran this function
% on a lowercase syllable (such as 'A', then 'a'), you will notice there is
% only one summary file *syl_A*.  If you run this function with 'a' then
% 'A' there is only one summary file *syl_a*.  This is an unavoidable bug
% in Windows operating system.  Recommendation: run a script to
% automatically relabel all of your uppercase syllables to a lowercase
% syllable that is not used to label any other syllable ('A' --> 'z').
%"Because the Windows operating system considers two files with the same name
%to be the same file (regardless of case), you cannot have two files with the 
%same name in the same folder. If you save MYFILE, and myfile.mat already 
%exists in the current folder, then MYFILE.MAT replaces myfile.mat without warning.
%If you save myfile, and MYFILE.mat already exists in the current folder, the
%contents of myfile.mat replace the contents of MYFILE.mat, but the name remains MYFILE.mat."
eval(sprintf(['save(''%s'',''P1_save_allfreqs'',''F1_save_allfreqs'',''P1_save'',''F1_save'',',...
    '''peak_pinterp'',''weighted_avg'',''amp_at_pitchquant'',''spect_entropy'',''wiener_entropy'',',...
    '''fname_arr'',''t_out'',''save_T1'','...
    '''spect_params'',''f_cutoff'',''f_cutoff_full'',''t_assay'',''t_assay_amp'',''birdname'',''syl_sum'')'],summary_fn));

end %end of function headphones_quantify_pitch

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Internal function to compute amplitude in a window around t_assay (ms
% since start of syllable)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [amp, window]=compute_summed_amp_window(smoothed,t_assay,Fs)
    window_half_width=8;   % in msec, the width of the window on EITHER side of t_assay
    window=t_assay+[-window_half_width window_half_width]*.001;% in seconds
    window_id=round(window*Fs);% id of time at start and end of window
    if window_id(1)<1;window_id(1)=1;end
    if window_id(2)>length(smoothed);window_id(2)=length(smoothed);end
    amp=mean(smoothed(window_id(1):window_id(2)));
    if isnan(amp);disp('NaN for amplitude, setting to zero');amp=0;end
    if length(amp)>1;error('length of amp variable > 1');end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Internal function to compute offset in seconds between sounds from 2
% channels by finding the lag where cross-correlation function peaks.
% Note: sm_sh could also be sm_mm.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function offset=compute_delay_from_xcorr(sm,sm_sh,Fs)
    [c,lags]=xcorr(sm_sh,sm);   % pos peak will be shifted lagging unshifted
    id_pos=find(lags>=0);   % only accept positive lags
    c=c(id_pos);
    lags=lags(id_pos);
    [~,peak_id]=max(c);
    offset=lags(peak_id)/Fs;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Internal function to parabolically interpolate 3 points and return
% x and y values at parabola's maximum.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xmax,ymax] = pinterp_internal(xs, ys)
    % Parabolically interpolates the peak given three points.
    % This is straightforward linear algebra.  You have 3 points (x1,y1),
    % (x2,y2), (x3,y3), and you want to find parabola that passes through these
    % points.  Parabola has equation y=a*x^2+b*x+c where a, b, c, are unknown.
    %Matrix X = x1^2 x1 1
    %           x2^2 x2 1
    %           x3^2 x3 1
    %Matrix Y = y1
    %           y2
    %           y3
    %Matrix ABC = a
    %             b
    %             c
    %System of 3 linear equations: X*ABC=Y --> inv(X)*X*ABC = inv(X)*Y --> ABC= inv(X)*Y
    %Matlab recommends x\ys instead of inv(x)*ys.
    x = [xs.^2, xs, [1;1;1]];
    abc =x\ys;

    % finds x and y where parabola peaks. xmax is the frequency that is saved to
    % file as the syllable's calculated pitch (peak_pinterp). ymax is the power
    % of this pitch.  The parabola's equation is y=a*x^2+b*x+c, to find
    % the maximum value you calculate where the derivative is 0.
    % y' = 2*a*x + b, find x where y'=0 --> 0 = 2*a*x + b --> xmax = -b/(2*a*x).
    % Once x is obtained, plug it in to get ymax = a*xmax^2 + b*xmax + c.
    xmax = -abc(2)/(2*abc(1));
    ymax = abc(1)*xmax^2+abc(2)*xmax+abc(3);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Internal function to compute pitch using parabolic interpolation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [peak_pinterp]=compute_scalar_pitch(spect_slice,F1_for_pitchquant,extra_input)
    if nargin==2;extra_input=0;extra_str=[];end

    % Compute a single scalar value of pitch based on spect_slice
    max_p_id=find(spect_slice==max(spect_slice));
    if sum(max_p_id==[1 length(spect_slice)])
        old_max_p_id=max_p_id;
        if max_p_id==1
            disp([extra_str 'Peak at LOWER extreme of range - looking for convex peak between limits'])
        elseif max_p_id==length(spect_slice)
            disp([extra_str 'Peak at UPPER extreme of range - looking for convex peak between limits'])
        end

        % Looks for a convex peak in between limits. How the expression works:
        % Assume for example spect_slice = [1 1 2 4 2 1 1] so there is a single
        % convex peak.  sign(diff(spect_slice)) = [0 1 1 -1 -1 0].  This identifies
        % where the derivative is positive or negative - where values are
        % increasing or decreasing. This tells you whether you are "walking uphill
        % or downhill from left to right".
        % diff(sign(diff(spect_slice))) = [1 0 -2 0 1]  This tells you where the
        % derivative switches from positive to negative.  If you are standing on
        % the peak, you were just walking uphill to get from the previous value
        % up to the peak (so that sign(diff(spect_slice)) was 1) and you can see
        % that you will walk downhill from peak to next value (sign(diff(spect_slice))
        % will be -1). Then at the peak, diff(sign(diff(spect_slice))) = -2.
        % Elsewhere, you will not be "walking up, peak, walking down", which means 
        % there is no peak and diff(sign(diff(spect_slice))) ~= -2
        % find(sign(diff(sign(diff(spect_slice))))<0) = 3. It asks at which index
        % is diff(sign(diff(spect_slice)))<0?  Remember that it is only < 0 if
        % there is a convex peak.  Here there is only one peak.  The last thing
        % is to add +1 (3+1 = 4 = index of the peak). The index 3 actually
        % points to the value right before the peak because of the way Matlab's 
        % diff function works.
        id_neg_2nd_der=find(sign(diff(sign(diff(spect_slice))))<0)+1;

        % Find largest magnitude peak.
        max_p_id=id_neg_2nd_der(find(spect_slice(id_neg_2nd_der)==max(spect_slice(id_neg_2nd_der))));
    end

    if isempty(max_p_id)% if NO convex peak found above, use old edge peak
        peak=F1_for_pitchquant(old_max_p_id);
        disp('NO CONVEX PITCH PEAK FOUND - USING ORIGINAL PEAK AT EDGE OF RANGE')
        peak_pinterp=peak;
    else % use convex peak
        x_pinterp=F1_for_pitchquant(max_p_id-1:max_p_id+1); % frequency at peak and either side of peak
        y_pinterp=spect_slice(max_p_id-1:max_p_id+1); % power at peak and either side of peak
        [peak_pinterp,~] = pinterp_internal(x_pinterp, y_pinterp); % parabolic interpolation
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Internal function to compute pitch using weighted average.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function weighted_avg=compute_scalar_pitch_weighted_avg(spect_slice,F1_for_pitchquant,extra_input)
    % finds peak pitch in the same way that evtaf_freq does
    % i.e. by taking a weighted average of frequencies around the maximum in
    % the power spectrum
    if nargin==2;extra_input=0;extra_str=[];end

    % Compute a single scalar value of pitch based on spect_slice
    max_p_id=find(spect_slice==max(spect_slice));
    if sum(max_p_id==[1 length(spect_slice)])
        old_max_p_id=max_p_id;
        if max_p_id==1
            disp([extra_str 'Peak at LOWER extreme of range - looking for convex peak between limits'])
        elseif max_p_id==length(spect_slice)
            disp([extra_str 'Peak at UPPER extreme of range - looking for convex peak between limits'])
        end

        % Looks for a convex peak in between limits
        %See comments in internal function compute_scalar_pitch
        id_neg_2nd_der=find(sign(diff(sign(diff(spect_slice))))<0)+1;   % this is length(spect_slice)-2
        max_p_id=id_neg_2nd_der(find(spect_slice(id_neg_2nd_der)==max(spect_slice(id_neg_2nd_der))));
    end

    if isempty(max_p_id) % if NO convex peak found above, use old edge peak
        weighted_avg=F1_for_pitchquant(old_max_p_id);
    else % use convex peak
        %In cases where the peak is near the beginning or end of the vectors,
        %then can't do weighted average for NBins_preferred bins on each side
        %because vectors stop. So do fewer bins.
        NBins_preferred=3;
        NBins=min([NBins_preferred max_p_id-1 length(F1_for_pitchquant)-max_p_id]);
        if NBins~=NBins_preferred;disp(['NBins = ' num2str(NBins) ' rather than preferred value of ' num2str(NBins_preferred)]);end
        peak_F1_vec=F1_for_pitchquant(max_p_id-NBins:max_p_id+NBins);
        peak_P1_vec=spect_slice(max_p_id-NBins:max_p_id+NBins);

        % Each frequency in the NBins_preferred bins on either side of the peak, 
        % including the peak bin, is weighted by its power.  Weighted sum is
        % the pitch.
        weighted_avg=sum(peak_F1_vec.*peak_P1_vec)/sum(peak_P1_vec);
    end
end
