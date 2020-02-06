function save_results(ts_tc, resultsFolder_path, filename)
% This function saves timeseries data into given txt filePath
%
% INPUTS:
%        ts_tc (1_by_1 timeseries.collection of timeseries): 
%        timeseries of collecion of timeseries data want to save
%        resultsFolder_path (N_by_1 char): path of results
%
% OUTPUTS:
%       Video object


% create list of time series
if isa(ts_tc,'timeseries')
    ts_list{1} = ts_tc;
else
    names = ts_tc.gettimeseriesnames;
    for n = 1:length(names)
        ts_list{n} = ts_tc.get(names{n});
    end
end 

header{1} = ['Time (', ts_tc.TimeInfo.Units, ')'] ;
data(:,1) = ts_tc.Time;  

for n = 1:length(ts_list)
    ts = ts_list{n};
    header{n+1} = [ ts.Name, ' (', ts.DataInfo.Units, ')'] ;
    data(:,n+1) = ts.Data;    
end


% Result folder path
if ~exist(resultsFolder_path, 'dir')
   mkdir(resultsFolder_path)
end
txt_filePath = [ resultsFolder_path, '/', filename, '.txt'];


% Save txt file
fid = fopen(txt_filePath,'w');
fprintf(fid,'%s\t', header{:});  
fprintf(fid,'\n');
dlmwrite(...
    txt_filePath, data, ...
    '-append', ...
    'delimiter', '\t', ...
    'precision', '%.4f');
fclose(fid);

end

