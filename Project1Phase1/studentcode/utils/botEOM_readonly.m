function sdot = botEOM_readonly(t, s, u, params)
% BOTEOM_READONLY Solve robot equation of motion
%   botEOM_readonly calculate the derivative of the state vector
%
% INPUTS:
% t      - 1 x 1, time
% 6      - 3 x 1, state vector = [x, y, theta]
% u      - 2 x 1, wheel speed output from controller [left right] (only used in simulation)
% params - struct, output from turtlebot3() and whatever parameters you want to pass in
%
% OUTPUTS:
% sdot   - 3 x 1, derivative of state vector s
%
% NOTE: You should not modify this function
% See Also: botEOM, turtlebot3

%************ EQUATIONS OF MOTION ************************
% Limit the wheel speeds due to actuator limits
u_limit = min([params.max_vel, params.max_omega*params.wheel_base]) / params.wheel_radius;

u_clamped = max(min(u, u_limit), -u_limit);

v = (u_clamped(1) + u_clamped(2)) * params.wheel_radius / 2;
w = (u_clamped(2) - u_clamped(1)) * params.wheel_radius / params.wheel_base;

% Compute velocities
theta = s(3);

xdot = v * cos(theta);
ydot = v * sin(theta);

thetadot = w;

% Assemble sdot
sdot = zeros(3,1);
sdot(1)  = xdot;
sdot(2)  = ydot;
sdot(3)  = thetadot;

end
