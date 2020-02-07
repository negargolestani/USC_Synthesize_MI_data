classdef FRAME < handle
    % This class represents a video frame
    properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        rgb                                                                 % RGB frame
    end
    
    methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % -----------------------------------------------------------------
        function [ self ] = FRAME(rgb)
            % Class constructor
            %       FRAME(rgb)
            %
            % INPUTS:
            %        rgb (N_by_M_by_3 double): rgb image
            %
            % OUTPUTS:
            %       FRAME object
            
            self.rgb = rgb;
        end
        % -----------------------------------------------------------------
        function [ objects ] = getobjects(self, varargin)
            % This function extracts objects from the frame
            %       self.getobjects(color_1, color_2, ...)
            %
            % INPUTS:
            %       color_i (1_by_1 char): color of targeted object
            %
            % OUTPUTS:
            %       objects (1_by_Ncolor cell array of struct array with             
            %       fields Centroid, BoundingBox, Orientation): detected
            %       objects
            
            medfilt_size = 5;
            binarize_threshold = 0.2;
            min_conctPixls = 300;
            
            objects = [];
            for n = 1:nargin-1
                diff_im = imsubtract( ...
                    self.rgb(:,:,strfind('rgb', varargin{n})), ...
                    rgb2gray(self.rgb));                                        % subtract given colored component from the grayscale image to extract the red components in the image.
                diff_im = medfilt2(diff_im, [medfilt_size medfilt_size] );      % use a median filter to filter out noise
                diff_im = imbinarize( diff_im, binarize_threshold );            % convert the resulting grayscale image into a binary image.
                diff_im = bwareaopen( diff_im, min_conctPixls );                % remove all those pixels less than 300px
                
                bw = logical( diff_im );                                        % label all the connected components in the image.
                
                objects{n} = regionprops(bw, ...
                    'Centroid', ...
                    'Orientation', ...
                    'BoundingBox') ;                                             % get property of regions
                % objects = [ objects(:)', objs(:)'];
            end
        end
        % -----------------------------------------------------------------
        function [ ] = show(self, varargin)
            % This function shows frame and object given as input (if any)
            %       self.show(objects_c1, objects_c2, ...)
            %
            % INPUTS:
            %       objects_ci (1_by_Nobjs cell array of objects): targeted
            %       object in the frame. each cell contains objects of a
            %       color 
            %
            % OUTPUTS:
            %       None
            
            imshow(self.rgb);
            if nargin>1
                objects = varargin{1};
                Ncolors = length(objects);
                LW = 4; L = 100; H = 30; color = jet(Ncolors);
                for n = 1:Ncolors
                    objs_color_n = objects{n};
                    for m = 1:length(objs_color_n)
                        drawline(objs_color_n(m), L, color(n,:), LW);
                        % drawrectangle(objs_color_n(m), L, H, color(n,:), LW);
                    end
                end
            end
        end
        % -----------------------------------------------------------------           
    end
end

