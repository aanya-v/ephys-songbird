function filescreen_fig_rhd

% No arguments, no return values
% 
% PURPOSE: Manually screen files for birdsong.
% This function lets you sort intan sound files into folders 
% based on their spectrograms.
% 
% INSTRUCTIONS:
% Hit "p" to play the audio recording
% Hit "s" to move to folder 'screened_song'
% Hit "c" to move to 'screened_call' folder
% Hit "k" to move to 'screened_spontaneous' folder
% Hit left arrow key to go back and correct a mistake.
% Hit left arrow key multiple times to go back farther.
% Hit "e" to exit and move the files you've flagged to move
% Files are deleted for real once you hit "k" or "s" on the last file.
% Close out the window to quit function. Files will not move.
%
% NOTES:
% * Assuming birdsong files created by Intan recording controller with .rhd
% extension

ext = 'rhd';
files=ls(['*.' ext]);
if strcmp(files,'')
    error(['No .' ext ' files in current directory: ' pwd]);
end
files=cellstr(files);
numfiles = size(files,1);
keepfiles = zeros(numfiles,1); %if keepfiles(f)==1, keep file f.  else delete file f.

figure
scrsz = get(0,'ScreenSize'); %below: x, y, width, height, relative to screen size
set(gcf,'Position',[scrsz(3)/100 scrsz(4)/4 scrsz(3)*(98/100),scrsz(4)/2]);

f=1;
while f<=numfiles %screen each file
    %display spectrogram of current file
    current_file = cell2mat(files(f));
    [fs,audio] = read_Intan_RHD2000_audio(current_file);
    FS=fs;
    Y=audio(1,:);
    spect_from_waveform_leila(audio,fs,1,[.9 32]); %sober lab function
%     clim([-21 -18])
%     clim([-10 -25])
    clim([-20 -5])
    colormap("turbo")
    set(gca,'ylim',[500 10000])
    set(gca,'fontsize',12,'fontweight','bold','color','k')
    plot_title = [strrep(files{f},'_','\_') '     ' num2str(f) ' of ' num2str(numfiles)];
    title(plot_title);
    
    if f==numfiles
        textbox = annotation('textbox',[0 0 1 1],'String',...
        'After this, files are deleted','fontsize',20,'fontweight','bold',...
        'HorizontalAlignment','center','VerticalAlignment','middle');
    end
    if f~=numfiles && exist('textbox','var')
        set(textbox,'Visible','off');
    end
    %wait until user hits 'k', 'd', left arrow key, 'p', 'e', 'c', or 's'
    c='';
    while ~strcmp(c,'k') && ~strcmp(c,'s')  && ~strcmp(c,'leftarrow') && ~strcmp(c,'p') && ~strcmp(c,'e') && ~strcmp(c,'c') %&& ~strcmp(c,'d')
        b=0;
        while b==0, b = waitforbuttonpress;  end
        c = lower(get(gcf, 'CurrentKey'));
    end
    if strcmp(c,'k') %spontaneous file (will move to 'screened_spont')
        title(plot_title,'color',[0.8 0 0]);
        keepfiles(f)=1;
    elseif strcmp(c,'s') %Mark as song file  (will move to 'screened_song')
        title(plot_title,'color',[0 .6 0]);
        keepfiles(f)=0;
%     elseif strcmp(c,'d') %delete file
%         title(plot_title,'color',[.8 0 0]);
%         keepfiles(f)=0;
    elseif strcmp(c,'c') %Mark as a call file
        title(plot_title,'color',[0 0.6 0]);
        keepfiles(f)=-1;

    %Allows user to correct mistakes on earlier files, such as pressing
    %'delete' when it should be 'keep'
    elseif strcmp(c,'leftarrow') %go to previous file.
        f = f-2; %this only works if outer loop is a WHILE loop
        if f==-1, f=0; end
        cla
    elseif strcmp(c,'p') % play song file
        soundsc(Y,FS)
        continue
    elseif strcmp(c,'e') % exit and delete files flagged for deleting
        keepfiles(f:end)=2; %keep the remaining files in the folder
        f=numfiles;
    end
    pause(.01)
    f = f+1;
end
close

%if blablabla.cbin was flagged for deleting, delete it, along with
% blablabla.tmp, blablabla.rec, and so on.

if exist('sreened_song','dir') == 0
    mkdir('screened_song')
end

if exist('calls','dir') == 0
    mkdir('calls')
end

if exist('screened_spont','dir') == 0
    mkdir('screened_spont')
end
% delfiles = files(keepfiles==0);
callfiles = files(keepfiles==-1);
songfiles = files(keepfiles==0);
spontfiles = files(keepfiles==1);

% if ~isempty(keepfiles)
%     disp('Keeping the following files.')
%     disp(char(keepfiles))
%     movefiles(keepfiles,'keep_files')
% end
if ~isempty(callfiles)
    disp('Move the following call files.')
    disp(char(callfiles))
    movefiles(callfiles,'calls')
end
if ~isempty(songfiles)
    disp('Move the following to screened_song folder.')
    disp(char(songfiles))
    movefiles(songfiles,'screened_song')
end
if ~isempty(spontfiles)
    disp('Move the following to screened_spont folder.')
    disp(char(spontfiles))
    movefiles(spontfiles,'screened_spont')
end
% if ~isempty(delfiles)
%    delfiles1 = strrep(delfiles, ext, '*');
%    cellfun(@delete, delfiles1);
% end
%    cd ..
%    cd mat\
%    delfiles2 = strrep(delfiles, ext, 'mat');
%    for j=1:length(delfiles2)
%        delfiles2{j}=['songdet1_' delfiles2{j}];
%    end
%    cellfun(@delete, delfiles2);
%    cd ..
%    cd gif\
%    delfiles3 = strrep(delfiles, ext, 'gif');
%    cellfun(@delete, delfiles3);
%    cd ..
%    cd wav\
end


function movefiles( cell_arr_of_filenames, destination, chunkSize )
% Move CHUNKSIZE files at a time. (system() commands can only be a certain character length)
if nargin<3, chunkSize = 100; end
n = length(cell_arr_of_filenames);
chunks = unique([(1:chunkSize:n) n]);
for i = 1:length(chunks)-1
    inds = chunks(i):chunks(i+1);
    fileListStr = sprintf(' "%s"', cell_arr_of_filenames{inds});
    cmd_str = sprintf('for %%f in (%s) do move %%f "%s"', fileListStr, destination);
    fprintf('Moving %d to %d of %d files to "%s" ...', inds([1 end]), n, destination)
    [status, result] = system(cmd_str);
    if status==0
        fprintf(' done!\n')
    else
        error(result)
    end
end
end