function [ synthMI ] = synthesize_MIdata(MIsys, motions)  
% This function synthesize the MI data for given motion of coils
%
%   INPUTS:
%           MIsys (1_by_1 MIsystem obj):  MIsystem object 
%           motions (1_by_Ncoils cell array of timeseries collection):
%           contain Location of coils' center and Normal of the coil's
%           surface 
%
%   OUTPUTS:
%          synthMI (1_by_Ntx cell array of timeseries): contains synthesize
%          MI data generated correspond to each TX coil (|S21(dB)|)

Ntx = min( length(MIsys.TX), length(motions)-1);
Nsamples = motions{1}.Length;

for t = 1:Nsamples
    RX_motion = motions{1};
    Crx = getdatasamples(RX_motion.Location, t);
    nrx = getdatasamples(RX_motion.Normal, t);
    if any([isnan(Crx),isnan(nrx)])
        S21_dB(t,:) = nan;
        continue
    else
        MIsys.RX.move (Crx, nrx);                                            % Move RX
    end
    
    for n = 1:Ntx
        TX_motion = motions{n+1};
        Ctx = getdatasamples(TX_motion.Location, t);
        ntx = getdatasamples(TX_motion.Normal, t);
        if any([isnan(Ctx),isnan(ntx)])
            S21_dB(t,:) = nan;
            continue
        else
            MIsys.TX(n).move (Ctx, ntx);                                    % Move TXs
        end
    end
    [ ~, S21_dB(t,:), ~, ~ ] = MIsys.SdB( );                                % Calculate S21 (dB) of TXs/Rx
end


% collection of timeseries object
time = motions{1}.Time;
synthMI = cell(1,Ntx);
for n = 1:Ntx
    synthMI{n} = timeseries( ...
        abs(fillmissing( S21_dB(:,n), 'linear')), ...                            % Clean data: remove NaNs
        time, ...
        'Name', ['TX_',num2str(n)] );
    synthMI{n}.DataInfo.Units = ' |S_21 (dB)| ';
    
    % Plot
    % figure, plot(synthMI{n}, 'LineWidth', 2);
    % xlim([synthMI{n}.TimeInfo.Start, synthMI{n}.TimeInfo.End]);
end
synthMI = tscollection(synthMI, 'Name','Synthetic MIdata');

end
