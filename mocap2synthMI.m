close all; clear all; clc;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Synthesize MI Motion Data (MoCap) %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization
% database = 'BML'; fileType = '.ptd';
database = 'MHAD'; fileType = '.txt';
mocapFolder_path = ['./data/mocap/',database];
resultsFolder_path = ['./results/synth_MImotion/',database];
addpath('src');

filename_list = lsfiles(mocapFolder_path, fileType);                        % List of all files in data directory

%% Execute code for all given filenames in filename_list
for n = 1:length(filename_list)    
    filename = filename_list{n};
    mocapFile_path = [mocapFolder_path,'/',filename];

    try
        synthMI_tscol = synthmimotion(database, mocapFile_path);
        savetstc(synthMI_tscol, resultsFolder_path, filename);               % save data
%         print(gcf, [resultsFolder_path,'/',filename], '-dtiff', '-r350');   % Save plot
   
    catch ME
        fprintf( 'Unsuccessful attempt: %s \n i = %d \t filename = %s \n %s \n ', ME.message, n, filename);
        continue;  % Jump to next iteration
    end
    
end