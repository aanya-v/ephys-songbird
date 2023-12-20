
tmp=pwd;
cd('/Users/ssober/Dropbox/Claire_tang/Metric_Space_Analysis (1)/datafiles/Premotor_window_40_0/C=2')
clear;load pu26y2_c_1415.mat
for x=1:length(spiketrains{1})
   count(x)=length(spiketrains{1}{x}); 
end

figure(1);clf;
subplot(2,1,1);hold on
plot(count,'ko','markerfacecolor','k')
xlabel('Trial #')
ylabel('# spikes / 40 msec')

subplot(2,1,2);hold on
plot(count,'k-o','markerfacecolor','k')
xlabel('Trial #')
ylabel('# spikes / 40 msec')

cd(tmp)