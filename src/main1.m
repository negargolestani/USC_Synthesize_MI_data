%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Measured vs. Simulated MI Data %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc;

%% Initialization --> modify 
filename = '16';
dataFolder_path = '../data';
resultsFolder_path = '../results';

%% Read files
video = Video([dataFolder_path, '/', filename, '.MOV']);                    % load Video
vna = VNA([dataFolder_path, '/', filename, '.s2p']);                        % VNA (timeseries)

%% Video: get motions
motions = video.get_motions( ...                                            % Get motions: translation and rotation of coils in the video
    2, ...                                                                  % Number of markers
    'r', ...                                                                % Marker color
    'b', ...                                                                % calLabel color
    0.2 ...                                                                 % Length of calLabel
    );

dist = motions{1}.Location - motions{2}.Location;                           % Distance between coils
dist.Data = fillmissing( dist.Data, 'linear');                              % Clean data: remove NaNs

misalign = motions{1}.Normal - motions{2}.Normal;                           % Misalignment between coils
misalign.Data = fillmissing( asin(misalign.Data(:,3)) * 180/pi, 'linear');  % Clean data: remove NaNs

%% Synthesize MI-motion data
TRX = Transceiver( ...                                                      % Transceivers (coil)
    Coil( 5e-2, 1, 10) ...                                                  % Coil (radius,  number of turns, wire AWG)
    );
MIsys = MIsystem( ...                                                       % MI systen: same transceivers are used for RX and TXs
    13.56e6, ...                                                            % operating frequency
    TRX, ...                                                                % reciever (RX)
    TRX, ...                                                                % transmitter (TXs)
    length(motions)-1);                                                     % number of TXs  (TXs are located around each bone described before)

synthMI_tscol = synthesize_MIdata(MIsys, motions);                          % synthetic MIdata (tscollection)
synthMI = synthMI_tscol.get(synthMI_tscol.gettimeseriesnames{1});           % synthetic MIdata (timeseries)
synthMI.Name = synthMI_tscol.Name;
synced_synthMI = sync(synthMI, vna);                                        % synchronize synthMI with vna data

%% Plot
figure,  set(gcf, 'Units', 'Inches', 'Position', [2,2,14,9]);

subplot(3,1,1)
plot(dist.Time, dist.Data, 'LineWidth', 2);
xlim([0, 30]);
leg = legend({'X', 'Y', 'Z'});
set(leg, 'Location', 'southeast');
title('Distance'); ylabel('meter')

subplot(3,1,2)
plot(misalign.Time, misalign.Data, 'LineWidth', 2);
xlim([0, 30]);
title('Misalignment'); ylabel('degree')

subplot(3,1,3), hold on
for ts = [vna, synthMI, synced_synthMI]
    plot(ts.Time, abs(ts.Data), ...
        'LineWidth', 2, ...
        'DisplayName', ts.Name);    
end
xlim([0, 30]);
leg = legend('show'); set(leg, 'Location', 'southeast');
title('Simulation vs Measurement'); 
ylabel('|S21 (dB)|')

%% Save results
save_results(synthMI, resultsFolder_path, filename);
save_results(synced_synthMI, resultsFolder_path, [filename, '_synced']);
% print(gcf, [resultsFolder_path,'/',filename], '-dtiff', '-r350');                        
