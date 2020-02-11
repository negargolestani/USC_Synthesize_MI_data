function [ filename_list ] = lsfiles(folder_path, type, varargin)
% This function returns list of all files with given type in the given
% directory path
%
% INPUTs:
%       folder_path (N_by_1 char): directory path
%       type (N_by_1 char): targeted type (e.g. ".txt")
%       removeType (boolean): remove type in filename or not
%
% OUTPUTS:
%       filename_list (N_by_1 cell array): name of all files in the given
%       folder

currDir = pwd();
cd( folder_path );
file_list = dir(['*',type]);
cd( currDir );


removeType = true;
if nargin == 3, removeType = varargin{1}; end    
  
% Remove type from end of names
if removeType
    N = length(file_list);
    filename_list = cell(N,1);
    l = length(type);
    for n = 1:N
        filename_list{n} = file_list(n).name(1:end-l);
    end    
else
    filename_list = {file_list.name};
end

end

