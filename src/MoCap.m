classdef MoCap <handle
    % This class represents Motion Capture (MoCap) data
    
    properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        X, Y, Z
        fs 
    end
    
    methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % -----------------------------------------------------------------
        function self = MoCap(varargin)
        % Class constructor
        %       MoCap(database, filepath)
        %       MoCap(X,Y,Z,fs)        
        %
        % INPUTs:
        %       database (N_by_1 char): name of databases
        %                               'BML', 'MHAD'
        %       filepath (N_by_1 char): directory of file        
        %       X, Y, Z (Ntime_by_Nmarkers double): markers xyz location (m)
        %       fs (1_by_1 double): sampling frequency (Hz)
        %
        % OUTPUTS:
        %       MoCap object

        if nargin == 2
            database = varargin{1};
            filepath = varargin{2};  
            
            % get filetype of given database
            switch database
                case 'BML'
                    filetype = '.ptd';
                case 'MHAD'
                    filetype = '.txt';                    
            end            
            if ~strcmp( filepath(end-3:end), filetype)
                filepath = [filepath, filetype];
            end
            
            % Read data
            [ self.X, self.Y, self.Z, self.fs ] = feval( ...
                ['readMoCap','_',database], filepath);            
            
        elseif nargin == 4
            self.X = varargin{1};
            self.Y = varargin{2};
            self.Z = varargin{3};
            self.fs = varargin{4};
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
        function [ new_self ] = downsample(self, M)
            % This function re-arrange/merge idxs of markers  
            %
            % INPUTS:
            %        M (1_by_1 double): downsample rate 
            %
            % OUTPUTS:
            %       new_self (1_by_1 MoCap obj): containing modified data
            
            new_self = MoCap( ...
                self.X(1:M:end,:), ...
                self.Y(1:M:end,:), ...
                self.Z(1:M:end,:), ...
                self.fs/M );
        end
        % -----------------------------------------------------------------
        function [ new_self ] = reshape_markers(self, markersIDX)
            % This function re-arrange/merge idxs of markers  
            %
            % INPUTS:
            %        markersIDX (Nmarkers_by_1 cell array): containing
            %        new_idx corresponding to each new idx 
            %
            % OUTPUTS:
            %       new_self (1_by_1 MoCap obj): containing modified data
            
            for n = 1:length(markersIDX)
                x(:,n) = nanmean( self.X(:, markersIDX{n}), 2);
                y(:,n) = nanmean( self.Y(:, markersIDX{n}), 2);
                z(:,n) = nanmean( self.Z(:, markersIDX{n}), 2);
            end
            new_self = MoCap(x, y, z, self.fs);
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