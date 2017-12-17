function [waypts] = circle
% CIRCLE generates waypoints to move the robot in a circle
% 
% The circle should have radius 1m and should be centered at the point 
% (0, 1)m.  The robot should traverse the trajectory in a counter-clockwise
% fashion.  The robot should start at the origin.
% 


r = 1;
n = 25;
theta = 3*pi/2;
Q = (n*10) / 21;
for i=1:21
    waypts(1,i) = r*cos(theta);
    waypts(2,i) = (r*sin(theta)+1);
    theta = theta + (pi/10);
end

% 
% waypts = [0,1,2;0,1,1];


% %define circle constraints
% radius = 1;
% center = [0 1];
% 
% %Number of sample points
% N = 50;
% theta=0:2*pi/N:2*pi;
% %theta = 0:0.5:2*pi;
% m = radius * [cos(theta')+center(1) sin(theta')+center(2)];
% 
% waypts = m';
% 
% 
% waypts = [0, waypts(1,:); 0, waypts(2,:)];
end
