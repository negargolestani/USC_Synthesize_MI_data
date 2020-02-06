function [ lag ] = get_lag(data1, data2)
% This function calculates lag of data1 respect to data2 using
% cross correlation
%
% INPUTS:
%       data1, data2 (1_by_1 timeseries): data can have different times
%
% OUTPUTS:
%       lag: lag of data1 with respect to data2 (s)


% Resample both data to time of data1 
Tmin = max(data1.TimeInfo.Start, data2.TimeInfo.Start);
Tmax = min(data1.TimeInfo.End, data2.TimeInfo.End);
time = data1.Time;
time = time(Tmin<time & time<Tmax);

data1_ = resample(data1, time);
data2_ = resample(data2, time);

% Find lag
c = corrcoef(data1_.Data, data2_.Data);
corr = [ c(1,2) ];
lag_range = [ 0 ];

% Correlation for different lags
for d = 1:50
    c = corrcoef( data1_.Data(1+d:end), data2_.Data(1:end-d) );
    corr = [corr, c(1,2)];
    lag_range = [lag_range, time(d)];
    
    c = corrcoef( data2_.Data(1+d:end), data1_.Data(1:end-d) );
    corr = [corr, c(1,2)];
    lag_range = [lag_range, -time(d)];
end
[lag_range, sortIdx] = sort(lag_range);
corr = corr(sortIdx);
lag = lag_range(corr == max(corr));


synced_data1 = data1;
synced_data1.Time = synced_data1.Time - lag;
end