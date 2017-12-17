function sdot = botEOM(t, s, controlhandle, trajhandle, params)
% botEOM Wrapper function for solving robot equation of motion
% 	botEOM takes in time, state vector, controller, trajectory generator
% 	and parameters and output the derivative of the state vector, the
% 	actual calculation is done in botEOM_readonly.
%
% INPUTS:
% t             - 1 x 1, time
% s             - 3 x 1, state vector = [x, y, theta]
% controlhandle - function handle of your controller
% trajhandle    - function handle of your trajectory generator
% params        - struct, output from turtlebot3() and whatever parameters you want to pass in
%
% OUTPUTS:
% sdot          - 3 x 1, derivative of state vector s
%
% NOTE: You should not modify this function
% See Also: botEOM_readonly, turtlebot3

% convert state to quad stuct for control
qd = stateToQd(s);

% Get desired_state
waypts = trajhandle();

% get control outputs
u = controlhandle(qd, waypts, t, params);

% compute derivative
sdot = botEOM_readonly(t, s, u, params);

end
