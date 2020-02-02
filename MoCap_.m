classdef MoCap
    % This class represents Motion Capture (MoCap) data
    
    properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        markers                                                             % timeseries collection of markers data
    end
    
    methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % -----------------------------------------------------------------
        function self = MoCap(filepath) 
            % Class constructor
            %       MoCap(filepath)
            %
            % INPUTs:
            %       filepath (1_by_N char): directory of MoCap data
            %
            % OUTPUTS:
            %       MoCap object
            %
            % NOTE:
            %       markers (1_by_15 cell array of Ntime_by_3 matrix)
            %       containing time-series xyz location of markers around 
            %       the body during movement. rows represent time and 
            %       columns represents xyz. 15 markers are placed as shown.
            %  
            %   (R)              N1
            %                    .
            %           N6 ......N2..... N3
            %          .    .         .    .
            %         N7    .         .     N4
            %        .      .         .      .
            %       N8     N13...N9....N10     N5
            %         	   .       .
            %         	   .       .
            %               N14      N11
            %         	   .       .
            %         	   .       .
            %               N15     N12
            %
            %       * To use another dataset, modify this function
            
            if strcmp(filepath(end-2:end),'ptd')
                self.markers = readMoCap_BML(filepath);
                
            elseif strcmp(filepath(end-2:end),'txt')
                self.markers = readMoCap_MHAD(filepath);            
            end                
        end
        % -----------------------------------------------------------------
        function [ ] = play(self)
            % This function shows MoCap data in 3D
            %
            % INPUTS:
            %        None
            %
            % OUTPUTS:
            %       None            
                        
            for n = 1:self.markers.size(2)
                markerName = ['Marker_',num2str(n)];
                marker = self.markers.get(markerName);
                
                X(:,n) = marker.Data(:,1);
                Y(:,n) = marker.Data(:,2);
                Z(:,n) = marker.Data(:,3);
            end
            xlimRange = [min(X(:)),max(X(:))];
            ylimRange = [min(Y(:)),max(Y(:))];
            zlimRange = [min(Z(:)),max(Z(:))];
            
            % Show data
            figure;
            for t = 1:size(X,1)
                x = X(t,:);
                y = Y(t,:);
                z = Z(t,:);
                scatter3(x, y, z, '*')
                xlim(xlimRange), ylim(ylimRange), zlim(zlimRange), drawnow
            end
        end
        % -----------------------------------------------------------------
        function [ motions ] = get_motion(self, bones_desc)
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
            time = self.markers.Time;
            
            for n = 1:Nbones
                N1 = self.markers.get(...
                    ['Marker_',num2str(bones_desc(n,1))]).Data;                
                N2 =  self.markers.get(...
                    ['Marker_',num2str(bones_desc(n,2))]).Data;
                
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ markers ] = readMoCap_BML(filepath)
% This function reads data of BML library and returns 15 markers. 
%
% INPUT:
%       fileName_action (N_by_1 char): file directory

% OUTPUT:
%       markers (collection of timeseries): collection of 15 timeseries
%       contain XYZ location of each markers (m)

%   * BML Library *
%   "Ma, Y., Paterson, & Pollick, F.E. (2006). A motion-capture library for 
%   the study of identity, gender, and emotion perception from biological 
%   motion. Behavior Research Methods, Instruments, & Computers, 38,
%   134-141."   
%   URL: http://paco.psy.gla.ac.uk/index.php/res/download-data


Nmarker = 15;

file = importdata(filepath);
Nt = file(1);
Nd = 45;
data = file(2:end);
data = reshape(data(1:Nt*Nd),Nd,Nt)';
data = data * 0.0254;                                                       % convert data from inch to meter
data = fillmissing( data, 'linear');                                        % clean data: remove NaNs

fs = 60;
time = [0:Nt-1] / fs;   

markers = cell(1,Nmarker);
for n = 1:Nmarker
    markers{n} = timeseries( ...
        data(: , 3*n-2:3*n), time, ...
        'Name', ['Marker_', num2str(n)] ...
        );
end
markers = tscollection(markers);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ markers ] = readMoCap_MHAD(filepath)
% This function reads data of MHAD library and returns 15 markers. 
%
% INPUT:
%       fileName_action (N_by_1 char): directory of file

% OUTPUT:
%       markers (collection of timeseries): collection of 15 timeseries
%       contain XYZ location of each markers (m)

Nmarker = 15;
M = 8;                                                                      % downsample
markers_idxs = {...
    1:3, ...
    [6,10], ...
    12:13,...
    14:15,...
    16:19,...
    20:21,...
    22:23,...
    24:27,...
    [4,7,8,11],...
    28:30,...
    31:32,...
    33:35,...
    36:38,...
    39:40,...
    41:43
    };

file = load(filepath );                                                     % load mocap data
data = file(1:M:end,1:end-2);                                               % Downsample
data = data * 1e-3;                                                         % convert data from cm to meter
data = fillmissing( data, 'linear');                                        % clean data: remove NaNs
Nt = size(data,1);

fs = 480/M;
time = [0:Nt-1] / fs;
for t = 1:Nt
    dt = reshape(data(t,:), 3, []);    
    
    % NaN all zero values
    for idx = find(sum(dt==0)==3)
        dt(:,idx ) = [nan, nan, nan];
    end
    Y(t,:) = dt(1,:);
    Z(t,:) = dt(2,:);
    X(t,:) = dt(3,:);    
end
 
X = fillmissing(X,'linear');
Y = fillmissing(Y,'linear');
Z = fillmissing(Z,'linear');

% timeseries 
markers = cell(1,Nmarker);
for n = 1:Nmarker    
    Xn = nanmean( X(:,markers_idxs{n}), 2);
    Yn = nanmean( Y(:,markers_idxs{n}), 2);
    Zn = nanmean( Z(:,markers_idxs{n}), 2); 
   
    markers{n} = timeseries( ...
        [Xn, Yn, Zn], time, ...
        'Name', ['Marker_', num2str(n)] ...
        );    
end
markers = tscollection(markers);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
