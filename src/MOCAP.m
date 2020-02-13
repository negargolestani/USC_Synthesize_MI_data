classdef MOCAP <handle
    % This class represents Motion Capture (MOCAP) data
    
    properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        X, Y, Z
        fs 
    end
    
    methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % -----------------------------------------------------------------
        function self = MOCAP(database, filepath)
        % Class constructor
        %       MOCAP(database, filepath)
        %
        % INPUTs:
        %       database (N_by_1 char): name of databases
        %                               'BML', 'MHAD'
        %       filepath (N_by_1 char): directory of file      
        %
        % OUTPUTS:
        %       MOCAP object
        
        [ X, Y, Z, self.fs ] = feval(['read',database], filepath);          % Read data
        err = find( (X==0).* (Y==0).* (Z==0) );
        X(err) = NaN;  Y(err) = NaN;  Z(err) = NaN;
        for n = 1:size(X,2)
            self.X(:,n) = fillmissing(X(:,n),'linear');
            self.Y(:,n) = fillmissing(Y(:,n),'linear');
            self.Z(:,n) = fillmissing(Z(:,n),'linear');
        end        
        
        end
        % -----------------------------------------------------------------
        function [ ] = play(self)
            % This function shows MOCAP data in 3D
            %
            % INPUTS:
            %        None
            %
            % OUTPUTS:
            %       None
            
            xlimRange = [min(self.X(:)),max(self.X(:))];
            ylimRange = [min(self.Y(:)),max(self.Y(:))];
            zlimRange = [min(self.Z(:)),max(self.Z(:))];
            
            % Show data
            figure;
            for t = 1:size(self.X,1)
                x = self.X(t,:);
                y = self.Y(t,:);
                z = self.Z(t,:);
                scatter3(x, y, z, '*')
                xlim(xlimRange), ylim(ylimRange), zlim(zlimRange), drawnow
            end
        end
        % -----------------------------------------------------------------
        function [ ] = downsample(self, M)
            % This function re-arrange/merge idxs of markers  
            %
            % INPUTS:
            %        M (1_by_1 double): downsample rate 
            %
            % OUTPUTS:
            %       new_self (1_by_1 MOCAP obj): containing modified data
            
            self.X = self.X(1:M:end,:);
            self.Y = self.Y(1:M:end,:);
            self.Z = self.Z(1:M:end,:);
            self.fs = self.fs/M ;
        end
        % -----------------------------------------------------------------
        function [ ] = merge_markers(self, markersIDX)
            % This function re-arranges/merges idxs of markers  
            %
            % INPUTS:
            %        markersIDX (Nmarkers_by_1 cell array): containing
            %        new_idx corresponding to each new idx 
            %
            % OUTPUTS:
            %       new_self (1_by_1 MOCAP obj): containing modified data
            
            for n = 1:length(markersIDX)
                x(:,n) = fillmissing( ...
                    nanmean( self.X(:, markersIDX{n}), 2), 'linear');
                y(:,n) = fillmissing( ...
                    nanmean( self.Y(:, markersIDX{n}), 2), 'linear');
                z(:,n) = fillmissing( ...
                    nanmean( self.Z(:, markersIDX{n}), 2), 'linear');
            end
            self.X = x; self.Y = y; self.Z = z;
        end
        % -----------------------------------------------------------------
        function [ motions ] = getmotion(self, bones_desc)
            % This function returns movement of requestd bones
            %
            % Inputs:
            %       bones_desc (Nbones_by_2 double): markers describing
            %       bones such that each bone is defined by two markers at
            %       its ends. rows represents bones and columns represents
            %       two markers at its ends.
            %
            % OUTPUTS:
            %       motions (Nbones_by_1 cell array of timeseries object):
            %       xyz location of bones center
            %       aligns (Nbones_by_1 cell array of Ntime_by_3 doubles):
            %       normalized alignment of bones
            
            Nbones = length(bones_desc);
            motions = cell(1,Nbones);                                       % collection of timeseries objects
            Nt = size(self.X,1);
            time = [0:Nt-1] / self.fs;   
            
            for n = 1:Nbones
                n1 = bones_desc(n,1);
                N1 = [self.X(:,n1),self.Y(:,n1),self.Z(:,n1) ];
                n2 = bones_desc(n,2);
                N2 = [self.X(:,n2),self.Y(:,n2),self.Z(:,n2) ];

                c_bone = (N1+N2)./2;                                        % Calculate the center of bones: assuming the coils are in the middle of bones
                loc = timeseries( ...
                    fillmissing( c_bone, 'linear'), time, ...               % Clean data: remove NaNs
                    'Name', 'Location' ...                                  % Centers location
                    );
                loc.DataInfo.Units = 'meter';
                
                n_bone = N2-N1;
                n_bone = n_bone ./( sqrt(sum(n_bone'.^2))' * ones(1,3));    % Calculate the alignment of line passing through two markers
                align = timeseries( ...
                    fillmissing( n_bone, 'linear'), time, ...               % Clean data: remove NaNs
                    'Name', 'Normal' ...                                    % Coils' surface normal
                    );
                
                motions{n} = tscollection( ...
                    {loc;align}, 'Name',['Coil_',num2str(n)] );
            end
            
        end
        % -----------------------------------------------------------------
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ X, Y, Z, fs ] = readBML(filepath, varargin)
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


filetype = '.ptd';
if ~strcmp( filepath(end-3:end), filetype)
    filepath = [filepath, filetype];
end
           
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
fs = 60;

end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ X, Y, Z, fs ] = readMHAD(filepath)
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

filetype = '.txt';
if ~strcmp( filepath(end-3:end), filetype)
    filepath = [filepath, filetype];
end

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
    X(t,:) = dt(1,:);
    Z(t,:) = dt(2,:);
    Y(t,:) = dt(3,:);    
end
fs = 480;

end