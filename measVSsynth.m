close all; clear all; clc;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Measured vs. Simulated MI Data %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization
videoFolder_path = './data/video';
vnaFolder_path = './data/vna';
resultsFolder_path = './results/meas_vs_synth';
addpath('src');

filename_list = intersect( ...
    lsfiles(videoFolder_path, '.MOV'), ...
    lsfiles(vnaFolder_path, '.s2p') ...
    );                                                                      % All filenames that exist for both video and VNA data

%% Execute code for all given filenames in filename_list
for n = 1:length(filename_list)
    filename = filename_list{n};
    
    [synced_synthMI, ~, ~, ~, ~] = measvssynth( ...
        videoFolder_path, vnaFolder_path, filename);
    
    print(gcf, [resultsFolder_path,'/',filename], '-dtiff', '-r350');       % Save plot
    savetstc(synced_synthMI, resultsFolder_path, filename);                 % save data
end