function [ vna ] = VNA( filepath )
% This function reads VNA data
%
% INPUTS:
%       filepath (N_by_1 char): directory of file
%
% OUTPUTS:
%       vna (1_by_1 timeseries): measured data by VNA (|S21(dB)|)

fid = fopen( filepath,'rt');
data = textscan(fid, '%f%f%f%f%f%f%f%f%f','HeaderLines',9);
fclose(fid);

time = data{1};
S21 = data{4}+1i*data{5};
S21_dB =  20*log10( abs(S21)) + 1i* (angle(S21)) ;

vna = timeseries( ...
    fillmissing( S21_dB, 'linear'), ...                                     % Clean data: remove NaNs
    time, ...
    'Name', 'VNA measurement');
vna.DataInfo.Units = ' S_21 (dB) ';

% Plot
% figure, plot(vna, 'LineWidth', 2);
% xlim([vna.TimeInfo.Start, vna.TimeInfo.End]);

end