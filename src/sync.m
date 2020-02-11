function [ synced_data1 ] = sync(data1, data2)
% This function synchronize data1 to data2 using get_lag function
%
% INPUTS:
%       data1, data2 (1_by_1 timeseries): data can have different times
%
% OUTPUTS:
%       synced_data1: synchronized data1 with data2 (timeseries)


lag = getlag(data1, data2);

time = data1.Time - lag;

synced_data1 = timeseries( ...
    data1.Data(time>0), ...
    time(time>0), ...
    'Name', ['Synchronized ', data1.Name] ...
    );
synced_data1.DataInfo.Units = data1.DataInfo.Units;

end