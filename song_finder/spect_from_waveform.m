function [S1,F1,T1,P1] =spect_from_waveform(Y,FS,plot,spect_params)
% [S1,F1,T1,P1] =spect_from_waveform(Y,FS,plot,spect_params)
% 
% SPECT_FROM_WAVEFORM  Make/plot spectrogram of smoothed waveform.  Give it 
% bin overlap and window size (ms) for fourier transform, default is
% [0 32].
% 
% ARGS:
% 
%     Y: raw waveform
%     FS: sampling rate Hz.
%     Get Y and FS by calling the function EVSOUNDIN.
%     
%     PLOT: 1 if you want to display a plot, 0 if you don't.
%     SPECT_PARAMS: [spect_overlap spect_win_dur], i.e.
%                   [percent_overlap window_size_in_milliseconds]
%                   default is [0 32]
% RETURNS:
% S1: spectrogram of the signal
% F1: frequencies at which spectrogram was computed
% T1: times at which the spectrogram was computed
% P1: power spectral density of each segment.  Natural log of P1 is
% plotted.

if nargin==2,plot=1;spect_params=[0 32];end
if nargin==3,spect_params=[0 32];end

spect_overlap = spect_params(1);  %percentage of overlap of specgram window
spect_win_dur = spect_params(2);


% BELOW GOOD FOR DISPLAY
%spect_overlap = 0.4;  %percentage of overlap of specgram window
%spect_win_dur=8;


% BELOW GOOD FOR ANALYSIS OF SYLLABLE 'C'
%spect_overlap = 0;  %percentage of overlap of specgram window
%spect_win_dur=32;
%disp(['spect program - Window duration ' num2str(spect_win_dur) ' msec, overlap = ' num2str(spect_overlap) ])
%spect_win_dur=16;
%first calculate nfft and noverlap
nfft=round(FS*spect_win_dur/1000);
% if FS>64000 & FS<65000 & nfft>1024 & nfft<1030
%     nfft=1024;
% end
% nfft = 2^nextpow2(nfft);
nfft=2^(round(log2(nfft))); % not up to next power of 2 - to nearest power of 2

%      disp('2x freq points')
% HACK IN *4 below

spect_win = hann(nfft);

noverlap = round(spect_overlap*length(spect_win)); %number of overlapping points

if 1%plot
    [S1,F1,T1,P1] =spectrogram(Y,spect_win,noverlap,nfft,FS);
else
    [S1,F1,T1,P1] =spectrogram_noplot(Y,spect_win,noverlap,nfft,FS);
end
%[nfft length(F1) F1(2:3)' diff(F1(2:3) )]

%[max(F1) length(F1)]
%length(Y)/FS

low_cutoff=500;
high_cutoff=12000;
%disp(['Only showing frequencies above ' num2str(low_cutoff) ' Hz'])
%disp(['Only showing frequencies below ' num2str(high_cutoff) ' Hz'])
id1=find(F1>=low_cutoff & F1<=high_cutoff);
F1=F1(id1);
P1=P1(id1,:);
S1=S1(id1,:);
mean(diff(F1));

if plot
    imagesc(T1,F1,log(P1)); set(gca,'YD','n');
    xlabel('Time (Seconds)'); ylabel('Frequency (Hz)'); 
end


