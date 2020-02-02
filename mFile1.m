%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Synthesize MI data - MoCap %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MoCap data Reference
%
%   * BML Library *
%   "Ma, Y., Paterson, & Pollick, F.E. (2006). A motion-capture library for 
%   the study of identity, gender, and emotion perception from biological 
%   motion. Behavior Research Methods, Instruments, & Computers, 38,
%   134-141."   
%   URL: http://paco.psy.gla.ac.uk/index.php/res/download-data

%   * MHAD Library *
%   " F. Ofli, R. Chaudhry, G. Kurillo, R. Vidal and R. Bajcsy. Berkeley 
%   MHAD: A Comprehensive Multimodal Human Action Database. In Proceedings
%   of the IEEE Workshop on Applications on Computer Vision (WACV), 2013."
%   URL: https://tele-immersion.citris-uc.org/berkeley_mhad

close all; clear all; clc;
%% Read files
mocap = MoCap('./data/ale_walk_nu_1_fin.ptd');                              % BML 
% mocap = MoCap('./data/moc_s01_a01_r01.txt');                                % MHAD

% mocap.play();                                                               % Play MoCap data

%% MoCap: get motions
motions = mocap.get_motion( ...                                             % Get motions: translation and rotation of coils around given bones 
    [ 9,2; 3,4; 4,5;10,11; 11,12; 6,7; 7,8; 3,14; 14,15] );                 % Description of human body skeleton bones: define each bone with two markers placed at its ends

%% MI systen: same transceivers are used for RX and TXs
TRX = Transceiver( ...                                                      % Transceivers (coil + matching)
    Coil( 5e-2, 1, 10), ...                                                 % Coil (radius,  number of turns, wire AWG)
    Inductor('S', 5.38e-6), ...                                             % Matching element (series inductor)
    Capacitor('P', 600e-12)...                                              % Matching element (parallel capacitor)
    );
MIsys = MIsystem( ...
    13.56e6, ...                                                            % Operating frequency
    TRX, ...                                                                % Reciever (RX)
    TRX, ...                                                                % Transmitter (TXs)
    length(motions)-1);                                                     % Number of TXs  (TXs are located around each bone described before)

% Synthesize MI-motion data
synthMI_tscol = synthesize_MIdata(MIsys, motions );

%% Plot
Ntx = length(motions)-1;
CM = jet(Ntx); FS = 10; LW = 2;

figure, hold on,  set(gcf, 'Units', 'Inches', 'Position', [4,0,10,10]);
for n = 1:Ntx
    synthMI = synthMI_tscol.get(synthMI_tscol.gettimeseriesnames{n});       % synthetic MIdata (timeseries)
    
    subplot(Ntx,1,n);
    plot( synthMI , ...
        'LineWidth', LW, ...
        'DisplayName', ['TX_{',num2str(n),'}'],...        
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