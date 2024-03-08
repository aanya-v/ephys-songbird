%% Sort rhd files into song-containing ('song') or spontaneous-only (spontaneous) folders. 
% Version by Leila, last edit: 7.13.22 
% Adapted from sort_rhd_lyndie
% No arguments, no return values
% 
% PURPOSE: Extracts and screens audio recordings for birdsong.
% This function lets you decide which sound files to keep and delete, 
% based on their spectrograms.
% 


function autosort_rhd
run_every_t=120;

tic
first=1;
while 1
    if first | toc>run_every_t
run_func
        pause(.1)
        tic
        first=0;
    end
end

function run_func

% make directories
if ~exist('spontaneous','dir')
    mkdir('spontaneous');
end

if ~exist('song','dir')
    mkdir('song');
end

% get all rhds except last one (which might still be open)
rhd=dir('*.rhd');
num_rhd = length(rhd);

% make sure there's not just one rhd (which would cause script to crash
% when it tries to get rid of the most recent rhd from the list
if num_rhd > 1
    rhd=rhd(1:end-1);
    % write list of rhds to text file 'batchfoo'
    fileID = fopen('batchfoo','w');
    for i=1:length(rhd)
        % '%s\n' = 'string newline'
        fprintf(fileID,'%s\n', rhd(i).name);
    end
    fclose(fileID)
else
    return
end

% Changed the last two variables here
cleandir_spect_SMART_intan('batchfoo',.020,1500,8,8); % first argument is the threshold - if sorting is not working properly try changing this value first

 
movebat_v2('batchfoo.dcrd','spontaneous')
%delbat_v2('batchfoo.dcrd')



movebat_v2('batchfoo.keep','song')
