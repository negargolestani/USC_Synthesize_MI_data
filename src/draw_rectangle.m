function [ ] = draw_rectangle(object, L, H, color, LW)
% This function plots rotated rectangles
%
% INPUTS:
%       object (1_by_1 struct): targeted objects that want to show on frame
%       L (1_by_1 double): length of rectangle (pixles)
%       H (1_by_1 double): height of rectangle (pixles)
%       color (Matlab color): color of lines
%       LW (1_by_1 double): width of lines
%
% OUTPUTS:
%       None

center = object.Centroid;
theta = object.Orientation * pi/180;

R = ([ cos(theta), -sin(theta); -sin(theta), -cos(theta) ]);

X = ([-L/2, L/2, L/2, -L/2]) ;
Y = ([-H/2, -H/2, H/2, H/2]) ;

f = [1, 2, 3, 4];
v = ( R * [X; Y] )';
for i = 1:2
    v(:,i) = v(:,i) + center(i);
end

patch('Faces', f,...
    'Vertices', v,...
    'FaceColor','none',...
    'Edgecolor',color,...
    'Linewidth',LW);
axis equal;

end