function [ synced_synthMI, synthMI, measuredMI, dist, misalign ] = measvssynth( ...
    videoFolder_path, vnaFolder_path, filename)
% This function synthesize the MI data for given motion of coils and
% synchronize with its corresponding vna data.
%
%   INPUTS:
%           videoFolder_path (N_by_1 char): path directory of folder
%                                           containing camera videos
%           videoFolder_path (N_by_1 char): path directory of folder
%                                           containing vna data
%           filename (N_by_1 char): name of file (must be same for both vna
%                                            and video files)
%
%   OUTPUTS:
%          synced_synthMI (1_by_1 timeseries): synchronize synthetic MI
%          data with measure data (S21dB)
%          synthMI (1_by_1 timeseries): synthetic MI data (S21dB)
%          measuredMI (1_by_1 timeseries): measured MI data (S21dB)
%          dist (1_by_1 timeseries): relative distance between coils
%          misalign (1_by_1 timeseries): relative misalignment between coils


% Video -> motions
video = VIDEO([videoFolder_path, '/', filename, '.MOV']);                   % load Video
motions = video.getmotion( ...                                              % Get motions: translation and rotation of coils in the video
    2, ...                                                                  % Number of markers
    'r', ...                                                                % Marker color
    'b', ...                                                                % calLabel color
    0.2 ...                                                                 % Length of calLabel
    );

% MI system
TRX = TRANSCEIVER( ...                                                      % Transceivers (coil)
    COIL( 5e-2, 1, 10) ...                                                  % Coil (radius,  number of turns, wire AWG)
    );
MIsys = MISYSTEM( ...                                                       % MI systen: same transceivers are used for RX and TXs
    13.56e6, ...                                                            % operating frequency
    TRX, ...                                                                % reciever (RX)
    TRX, ...                                                                % transmitter (TXs)
    length(motions)-1);                                                     % number of TXs  (TXs are located around each bone described before)

% Synthetic MI data
synthMI_tscol = MIsys.synthesizedata(motions);                                  % synthetic MIdata (tscollection)
synthMI = synthMI_tscol.get(synthMI_tscol.gettimeseriesnames{1});           % synthetic MIdata (timeseries)
synthMI.Name = synthMI_tscol.Name;

% Measured MI data (VNA)
vnaData = readvna([vnaFolder_path, '/', filename, '.s2p']);                 % VNA (timeseries)
measuredMI = vnaData.S21;
measuredMI.Name = 'Measured MI data';

% Synchronize synthetic data with measurements
synced_synthMI = sync(synthMI, measuredMI);                                 % synchronize synthMI with vna data


% Relative movement
dist = motions{1}.Location - motions{2}.Location;                           % Distance between coils
dist.Data = fillmissing( dist.Data, 'linear');                              % Clean data: remove NaNs
dist.Name = 'Distance';

misalign = motions{1}.Normal - motions{2}.Normal;                           % Misalignment between coils
misalign.Data = fillmissing( asin(misalign.Data(:,3)) * 180/pi, 'linear');  % Clean data: remove NaNs
misalign.Name = 'Misalignment';
misalign.DataInfo.Units = 'degree';


% Plot results
plot_(synced_synthMI, synthMI, measuredMI, dist, misalign);

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_(synced_synthMI, synthMI, measuredMI, dist, misalign)
LW = 2;

figure('Name','Simulation vs Measurement','NumberTitle','off');
set(gcf, 'Units', 'Inches', 'Position', [2,2,14,9]);

subplot(3,1,1)
plot(dist.Time, dist.Data, 'LineWidth', LW);
xlim([0, 30]);
leg = legend({'X', 'Y', 'Z'});
set(leg, 'Location', 'southeast');
title(dist.Name);
ylabel(dist.DataInfo.Units)

subplot(3,1,2)
plot(misalign.Time, misalign.Data, 'LineWidth', LW);
xlim([0, 30]);
title(misalign.Name);
ylabel(misalign.DataInfo.Units)

subplot(3,1,3), hold on
for ts = [measuredMI, synthMI, synced_synthMI]
    plot(ts.Time, abs(ts.Data), ...
        'LineWidth', LW, ...
        'DisplayName', ts.Name);
end
xlim([0, 30]);
leg = legend('show'); set(leg, 'Location', 'southeast');
title('Simulation vs Measurement');
ylabel('|S21|  (dB)' )
end