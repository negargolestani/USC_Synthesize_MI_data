function [ synthMI_tscol ] = synthmimotion( database, mocapFile_path)
% This function synthesize the MI motion data for given motion of coils and
% synchronize with its corresponding vna data.
%
%   INPUTS:
%           mocapFile_path (N_by_1 char): directory of file      
%           database (N_by_1 char): name of database 'BML', 'MHAD'
%
%   OUTPUTS:
%          synthMI_tscol (1_by_1 collection of timeseries): generated
%          synthetic MI motion data corresponding to given file for eight
%          TX around human body  (S21dB)         
% ---------------------------------------------------
%       standard 15-marker
% ---------------------------------------------------
%   (R)              N1
%                    .
%           N6 ......N2..... N3
%          .    .         .    .
%         N7    .         .     N4
%        .      .         .      .
%       N8     N13...N9....N10     N5
%         	   .       .
%         	   .       .
%               N14      N11
%         	   .       .
%         	   .       .
%               N15     N12
% ---------------------------------------------------



% MoCap -> motions
mocap = MOCAP(database, mocapFile_path);

% convert to  standard 15-marker
if strcmp(database, 'MHAD')
    Nodes_idx = { 1:3, [6,10], 12:13, 14:15, 16:19, 20:21, 22:23, ... 
        24:27, [4,7,8,11], 28:30,  31:32, 33:35, 36:38, 39:40, 41:43 };
    mocap.merge_markers(Nodes_idx)
end

% mocap.play();                                                             % Play MoCap data
motions = mocap.getmotion( ...                                              % Get motions: translation and rotation of coils around given bones
    [ 9,2; 3,4; 4,5;10,11; 11,12; 6,7; 7,8; 3,14; 14,15] );                 % Description of human body skeleton bones: define each bone with two markers placed at its ends

% MI system
TRX = TRANSCEIVER( ...                                                      % Transceivers (coil + matching)
    COIL( 5e-2, 1, 10), ...                                                 % Coil (radius,  number of turns, wire AWG)
    INDUCTOR('S', 5.38e-6), ...                                             % Matching element (series inductor)
    CAPACITOR('P', 600e-12)...                                              % Matching element (parallel capacitor)
    );
MIsys = MISYSTEM( ...
    13.56e6, ...                                                            % Operating frequency
    TRX, ...                                                                % Reciever (RX)
    TRX, ...                                                                % Transmitter (TXs)
    length(motions)-1);                                                     % Number of TXs  (TXs are located around each bone described before)

% Synthetic MI data
synthMI_tscol = MIsys.synthesizedata(motions);                              % synthetic MIdata (tscollection)

% Plot results
plot_(synthMI_tscol)

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_(synthMI_tscol)

names = synthMI_tscol.gettimeseriesnames;
Ntx = length(names);
CM = jet(Ntx); 
FS = 10; 
LW = 2;

figure('Name','Synthetic MI Motion','NumberTitle','off');
set(gcf, 'Units', 'Inches', 'Position', [4,0,10,10]);
hold on,  
for n = 1:Ntx
    synthMI = synthMI_tscol.get(names{n});                                  % synthetic MIdata (timeseries)
    
    subplot(Ntx,1,n);
    plot(synthMI.Time, abs(synthMI.Data) , ...
        'LineWidth', LW, ...
        'DisplayName', names{n},...
        'Color', CM(n,:) ...
        );
    
    leg = legend(gca,'show');
    set(leg, 'FontSize',FS, 'loc','eastoutside');
    
    ytickformat('%.1f');
    xlim([synthMI.TimeInfo.Start, synthMI.TimeInfo.End])
    title('');
    ylabel('');
    if n<Ntx
        set(gca,'XTick',[]);
        xlabel('');
    end
end

end