function [ robot ] = robot_pos( pos, rot, r, L, R, n )
%ROBOT_POS Calculates coordinates of robot's position in world frame
% pos       2x1 position vector [x; y];
% rot       2x2 body-to-world rotation matrix
% r         1x1 wheel radius
% L         1x1 wheel base
% R         1x1 robot radius
% n         1x1 number of points for circle

if nargin < 6
    n = 20;
end

wHb   = [rot pos(:); 0 0 1]; % homogeneous transformation from body to world

robotBodyFrame  = [r L/2 1; -r L/2 1; r -L/2 1; -r -L/2 1]';

th = linspace(0, 2*pi, n);
robotBodyFrame  = [robotBodyFrame, [R*cos(th); R*sin(th); ones(1,n)]];

robotWorldFrame = wHb * robotBodyFrame;
robot           = robotWorldFrame(1:2, :);

end
