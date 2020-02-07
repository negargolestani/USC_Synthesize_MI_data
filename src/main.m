close all; clear all; clc;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Measured vs. Simulated MI Data %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Directories 
videoFolder_path = '../data/video';
vnaFolder_path = '../data/vna';
resultsFolder_path = '../results/meas_vs_synth';
if ~exist(resultsFolder_path, 'dir'), mkdir(resultsFolder_path); end

filename = '16';

[synced_synthMI, ~, ~, ~, ~] = measvssynth(videoFolder_path, vnaFolder_path, filename);
print(gcf, [resultsFolder_path,'/',filename], '-dtiff', '-r350');           % Save plot
savetstc(synced_synthMI, resultsFolder_path, filename);                     % save data


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Synthesize MI Motion Data (MoCap) %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mocapFolder_path = '../data/mocap';
resultsFolder_path = '../results/synth_MImotion';
if ~exist(resultsFolder_path, 'dir'), mkdir(resultsFolder_path); end

% database = 'BML';  filename = 'ale_walk_nu_1_fin';                          
database = 'MHAD'; filename = 'moc_s01_a01_r01';

synthMI_tscol = synthmimotion( mocapFolder_path, filename, database);
print(gcf, [resultsFolder_path,'/',filename], '-dtiff', '-r350');           % Save plot
savetstc(synthMI_tscol, resultsFolder_path, filename);                      % save data
