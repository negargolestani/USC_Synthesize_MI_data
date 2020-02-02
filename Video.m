classdef Video
    
    properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        vid                                                                 % VideoReader
    end
    
    methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % -----------------------------------------------------------------
        function [ self ] = Video(filepath)
            % Class constructor
            %       Video(filepath)
            %
            % INPUTS:
            %        filepath (N_by_1 char): directory of video file
            %
            % OUTPUTS:
            %       Video object
            
            self.vid = VideoReader(filepath);
        end
        % -----------------------------------------------------------------
        function [ frame ] = next(self)
            % This function gets next frame of video
            %
            % INPUTS:
            %       None
            %
            % OUTPUTS:
            %       frame (1_by_1 Frame obj): next frame of video
            
            if hasFrame(self.vid)
                frame = readFrame(self.vid);
                frame = Frame(frame);
            else
                frame = false;
            end
        end
        % -----------------------------------------------------------------
        function [  ] = play(self, duration, varargin)
            % This function plays video and shows object with targeted
            % colod (if any given)
            %
            % INPUTS:
            %       duration (1_by_1 double): duration of play time (s)
            %       color (1_by_1 char): color of targeted objects
            %
            % OUTPUTS:
            %       None
            
            while hasFrame(self.vid)
                frame = self.next();
                frame.show( varargin{:}  )
                pause(0.01);
                
                if self.vid.CurrentTime > duration
                    break
                end
                
            end
        end
        % -----------------------------------------------------------------
        function [ pix2m_scale ] = calibrate(self, calLabel_color, calLabel_length_m)
            % This function converts pixles to meter using a calLabel
            % (calibration Label)
            % 
            % INPUTS:
            %       calLabel_color (1_by_1 char): color of calLabel 
            %       calLabel_length_m (1_by_1 double): length of calLabel (m)
            %
            % OUTPUTS:
            %       pix2m_scale (1_by_1 double): scale to convert pixles to
            %       meter

            frame  = self.next();
            calLabel = frame.get_objects( calLabel_color );
            calLabel = calLabel{1};
            
            theta = sign(calLabel.Orientation) * (90 - abs(calLabel.Orientation) );
            new_frame = Frame( imrotate(frame.rgb, theta) );
            new_calLabel = new_frame.get_objects( calLabel_color );
            new_calLabel = new_calLabel{1};
            calLabel_length = new_calLabel.BoundingBox(end);
            
            pix2m_scale = calLabel_length_m / calLabel_length;
        end
        % -----------------------------------------------------------------
        function [ motions ] = get_motions(self, ...
                Nmarker, marker_color, calLabel_color, calLabel_length_m)
            % This dunction returns motions of coils in the video
            %
            % INPUTS:
            %       Nmarkers (1_by_1 double): number of markers (coils)
            %       marker_color (1_by_1 char): color of markers
            %       calLabel_color (1_by_1 char): color of calLabel 
            %       calLabel_length_m (1_by_1 double): length of calLabel (m)
            %
            % OUTPUTS:
            %       motions (1_by_2 cell array of timeseries collections):
            %       contains the motions (Location and Normal) of coils
            
            
            pix2m_scale = self.calibrate(calLabel_color,calLabel_length_m); % Video calibration (using calLabel)
            locs = cell(1,Nmarker);
            aligns = cell(1,Nmarker);
            
            while hasFrame(self.vid)
                frame = self.next();
                coils = frame.get_objects( marker_color);
                coils = coils{1};
                xy = [coils.Centroid]; y = xy(2:2:end);
                [~, sortedidx] = sort(y);
                coils = coils(sortedidx);                                   % Sort coils based on their height
                
                if length(coils) == Nmarker
                    for n = 1:Nmarker
                        C = coils(n).Centroid * pix2m_scale;
                        locs{n}(end+1,:) = [0, C(1), C(2)];
                        alpha = ( coils(n).Orientation +90 )* pi/180;
                        aligns{n}(end+1,:) = [0, cos(alpha), sin(alpha)];
                    end
                else
                    for n = 1:Nmarker
                        locs{n}(end+1,:) = [nan, nan, nan];
                        aligns{n}(end+1,:) = [nan, nan, nan];
                    end
                end
            end
            
            motions = cell(1,Nmarker);                                      % collection of timeseries objects
            time = [0:length(locs{1})-1] / self.vid.FrameRate;
            for n = 1:Nmarker
                loc = timeseries( ...
                    fillmissing( locs{n}, 'linear'), time, ...              % Clean data: remove NaNs
                    'Name', 'Location' ...                                  % Centers location
                    );
                loc.DataInfo.Units = 'meter';
                align = timeseries( ...
                    fillmissing( aligns{n}, 'linear'), time, ...            % Clean data: remove NaNs
                    'Name', 'Normal' ...                                    % Coils' surface normal
                    );
                motions{n} = tscollection( ...
                    {loc;align}, 'Name',['Coil_',num2str(n)] );
                
                % Plot
                % figure,
                % subplot(2,1,1),plot(loc, 'LineWidth', 2);
                % xlim([loc.TimeInfo.Start, loc.TimeInfo.End]);
                % legend({'X', 'Y', 'Z'});
                % subplot(2,1,2), plot(align, 'LineWidth', 2);
                % xlim([align.TimeInfo.Start, align.TimeInfo.End]);
                % legend({'X', 'Y', 'Z'});
            end
            
        end
        % -----------------------------------------------------------------
        
    end
end