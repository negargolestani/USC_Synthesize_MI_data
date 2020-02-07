function [ vnaData ] = readvna( filepath )
% This function reads VNA data
%
% INPUTS:
%       filepath (N_by_1 char): directory of file
%
% OUTPUTS:
%       vnaData (1_by_1 collection of timeseries): S-parameters measured ,,,
%       by VNA 

fid = fopen( filepath,'rt');
data = textscan(fid, '%f%f%f%f%f%f%f%f%f','HeaderLines',9);
fclose(fid);

time = data{1};
names = {'S11','S21','S12','S22'};

for n = 1:4
    s = data{2*n}+1i*data{2*n+1};
    sdB = 20*log10( abs(s)) + 1i*(angle(s) );    
    sdB_ts = timeseries( ...
        fillmissing( sdB, 'linear'), ...
        time, 'Name', names{n});
        sdB_ts.DataInfo.Units = 'dB';
    
    SdB_list{n} = sdB_ts;
end

vnaData = tscollection( SdB_list, 'Name', 'VNA measurement');                   

end