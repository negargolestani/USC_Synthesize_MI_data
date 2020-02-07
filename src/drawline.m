function [ ] = drawline(object, L, color, LW)
% This function plots rotated line
%
% INPUTS:
%       object (1_by_1 struct): targeted objects that want to show on frame
%       L (1_by_1 double): length of line (pixles)
%       color (Matlab color): color of lines
%       LW (1_by_1 double): width of lines
%
% OUTPUTS:
%       None

C = object.Centroid;
alpha = object.Orientation * pi/180;

dx = L * cos(alpha);
dy = L * sin(alpha);

X = C(1) + [dx, -dx];
Y = C(2) + [-dy, dy];
line( X, Y, 'Color',color, 'LineWidth', LW )

% plot(C(1),C(2), '-w+')
% a = text( C(1), C(2)-30, ...
%     strcat('X: ', num2str(round(C(1))),...
%     '    Y: ', num2str(round(C(2)))));
% set( a,'Color', color);
end
