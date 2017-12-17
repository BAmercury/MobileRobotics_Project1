function [u] = controller_one(qd, waypts, t, params)
% CONTROLLER robot controller
%
% INPUTS:
% qd     - qd structure, contains current state
% waypts - 2 x n, way points [x; y]
% t      - 1 x 1, time
% params - struct, output from turtlebot3() and whatever parameters you want to pass in
%
% OUTPUTS:
% u      - 2 x 1, wheel velocities [left, right]

% =================7== Your code goes here ===================

lookahead = 0.09;
%lookahead = 2.00;



[pt_min, dist_min, segment_index] = find_closest_point(qd.pos,waypts);
   

world_goal = [0;0;1];


controller = robotics.PurePursuit;

controller.LookaheadDistance = 5;
controller.MaxAngularVelocity = params.max_omega;
controller.DesiredLinearVelocity = params.max_vel;
controller.Waypoints = transpose(waypts);


[v, w] = step(controller,[qd.pos;0]);

c_input = [v;w];

%% Set wheel velocities given actuator limits and relationship R = v/w

control_t = [(1/params.wheel_radius), ((-params.wheel_base) / (2*params.wheel_radius)); (1/params.wheel_radius), ((params.wheel_base) / (2* params.wheel_radius))];

u = control_t * c_input;





            

        




% =================== Your code ends here ===================


end
