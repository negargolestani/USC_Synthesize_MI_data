function [ X, Y, Z, fs ] = readMoCap_MHAD(filepath)
% This function reads data of MHAD library and returns all markers. 
%
% INPUT:
%       fileName_action (N_by_1 char): directory of file

% OUTPUT:
%       X,Y,Z (Ntime_by_Nmarkers double): location of markers (m)
%       fs (1_by_1 double): sampling frequency (Hz)
%
%   * MHAD Library *
%   " F. Ofli, R. Chaudhry, G. Kurillo, R. Vidal and R. Bajcsy. Berkeley 
%   MHAD: A Comprehensive Multimodal Human Action Database. In Proceedings
%   of the IEEE Workshop on Applications on Computer Vision (WACV), 2013."
%   URL: https://tele-immersion.citris-uc.org/berkeley_mhad

file = load(filepath );                                                     % load mocap data
data = file(:,1:end-2);                                                     % Downsample
data = data * 1e-3;                                                         % convert data from cm to meter
data = fillmissing( data, 'linear');                                        % clean data: remove NaNs
Nt = size(data,1);

for t = 1:Nt
    dt = reshape(data(t,:), 3, []);    
    for idx = find(sum(dt==0)==3)
        dt(:,idx ) = [nan, nan, nan];                                       % NaN all zero values                         
    end
    Y(t,:) = dt(1,:);
    Z(t,:) = dt(2,:);
    X(t,:) = dt(3,:);    
end
X = fillmissing(X,'linear');
Y = fillmissing(Y,'linear');
Z = fillmissing(Z,'linear');
fs = 480;
end

