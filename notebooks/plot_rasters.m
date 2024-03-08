birdID= 'br177yw112';
neuronID = 'day68';
syl_or_seq = 'c';
neural_channel = 2;
pre_onset_time_s = 0.6;
filelist = 'batch_passed_files.txt';
compile_case_v2(birdID,neuronID,syl_or_seq,neural_channel,pre_onset_time_s,filelist)
%% 
motif_sumfile = 'br177yw112_c_day62_ch2_premotor_40ms_spiketimes_acoustics_2023-12-20.mat';
song_file = 'br177yw112_221127_195302_songbout1.mat';
align_syl_num_all = 1;
which_trials = 'all';
trial_order = 'chron';
burstDetector = 0;
plot_multiple_FRs(song_file, motif_sumfile,align_syl_num_all,which_trials,trial_order, burstDetector)