function params = turtlebot3()
% turtlebot3: physical parameters for the Turtlebot 3 Burger
%
% 2017 Philip Dames
%
% This function creates a struct with the basic parameters for the
% Turtlebot 3 Burger differential drive robot
%

params.mass         = 1.0;   % weight (in kg) with battery, sensors, ets
params.wheel_base   = 0.160; % wheel base m
params.wheel_radius = 0.033; % wheel radius m
params.radius       = 0.105; % robot radius m
params.max_vel      = 0.22;  % maximum translational speed m/s
params.max_omega    = 2.84;  % maximum rotational speed rad/s
params.max_wheel    = 6.88;

% You can add any fields you want in params
% for example you can add your controller parameters by
% params.k = 0, and they will be passed into controller.m

end
