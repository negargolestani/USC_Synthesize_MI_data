function [ X, Y, Z, fs ] = readMoCap_BML(filepath)
% This function reads data of BML library and returns 15 markers. 
%
% INPUT:
%       filepath (N_by_1 char): file directory

% OUTPUT:
%       X,Y,Z (Ntime_by_Nmarkers double): location of markers (m)
%       fs (1_by_1 double): sampling frequency (Hz)
%
%   * BML Library *
%   "Ma, Y., Paterson, & Pollick, F.E. (2006). A motion-capture library for 
%   the study of identity, gender, and emotion perception from biological 
%   motion. Behavior Research Methods, Instruments, & Computers, 38,
%   134-141."   
%   URL: http://paco.psy.gla.ac.uk/index.php/res/download-data


file = importdata(filepath);
Nt = file(1);
Nd = 45;
data = file(2:end);
data = reshape(data(1:Nt*Nd),Nd,Nt)';
data = data * 0.0254;                                                       % convert data from inch to meter
data = fillmissing( data, 'linear');                                        % clean data: remove NaNs

X = data(:,1:3:end);
Y = data(:,2:3:end);
Z = data(:,3:3:end);

X = fillmissing(X,'linear');
Y = fillmissing(Y,'linear');
Z = fillmissing(Z,'linear');
fs = 60;
end