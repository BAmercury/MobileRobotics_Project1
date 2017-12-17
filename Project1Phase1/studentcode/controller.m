function [u] = controller(qd, waypts, t, params)
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
lookahead = 0.1;



[pt_min, dist_min, segment_index] = find_closest_point(qd.pos,waypts);
   

world_goal = [0;0;1];


%% Finding the Goal Point



% Find the endpoint of the current segment
end_point = waypts(:,(segment_index+1));
base_point = waypts(:,segment_index);

% Does the endpoint of current segment end OUTSIDE of the lh circle?
dist2endpoint = qd.pos - end_point;
dist2endpoint = norm(dist2endpoint);


if (dist2endpoint > lookahead)
    % Segment end is OUTSIDE circle, so closest point (pt_min) must be
    % inside. Goal point must be inside the point. Find it using trig
    %world_goal = [pt_min; 1];
    
   
    dist2goal = real(sqrt((lookahead^2) - (dist_min^2)));

    seg_direction = (end_point - base_point) / norm(end_point - base_point);
    
    
    
    point2goal = dist2goal * seg_direction;
    world_goal = [point2goal(1)+pt_min(1); point2goal(2)+pt_min(2);1];

    %take unit vector of segment to get direction, multiply the unit vector
    %by dist2goal
else
    % Check next segments until we find the one where the endseg is outside
    % lh
    
    start_segment = segment_index+1;
    max_segment = length(waypts);
    
    for (i = start_segment:max_segment)
        %search for a segment endpoint that lies outside lh circle
        end_point = waypts(:,i);
        dist2endpoint = qd.pos - end_point;
        dist2endpoint = norm(dist2endpoint);
        
        if (dist2endpoint > lookahead)
            endpoint_segment = i;
            break;
        else
            endpoint_segment = i;
        end

        
    end
    
    %Now we have the index of the endpoint. Let's find the goal point
    %base_segment = endpoint_segment - 1;
    
    [pt_min, dist_min] = find_closest_point(qd.pos,waypts, endpoint_segment);
    
    % Find the endpoint of the current segment
    if (endpoint_segment == length(waypts))
        end_point = waypts(:,endpoint_segment);
        base_point = waypts(:,(endpoint_segment-1));
    else
        
        end_point = waypts(:,(endpoint_segment+1));
        base_point = waypts(:,endpoint_segment);
    end
    
    dist2goal = real(sqrt( (lookahead^2) - (dist_min^2)));
    seg_direction = (end_point - base_point) / norm(end_point - base_point);



    
    point2goal = dist2goal * seg_direction;
    world_goal = [point2goal(1)+pt_min(1); point2goal(2)+pt_min(2);1];
    

    


        
end


%% Transform goal point

% Rotation of Robot Frame from World frame (wRr)
rot = [cos(qd.theta), -sin(qd.theta); sin(qd.theta), cos(qd.theta)];
% World frame from Robot Frame (rRw)
%rot = inv(rot);
%rot_h = rot;
% wTr
trans = [qd.pos];
%H = [rot_h, trans;0,0,1];
R = transpose(rot);

% wHr => rHw
H = [R, (-R * trans); 0,0,1];
robot_goal = H * world_goal


%% Calculate radius of the arc connecting robot to goal point

l_squared = (lookahead)^2;
R = l_squared / (2 * robot_goal(2));


%% Set wheel velocities given actuator limits and relationship R = v/w

control_t = [(1/params.wheel_radius), ((-params.wheel_base) / (2*params.wheel_radius)); (1/params.wheel_radius), ((params.wheel_base) / (2* params.wheel_radius))];


w = 1;



v = w * R;


%w = (v / R);



   
c_input = [v;w];





    
uint = control_t * c_input;


if (R == inf) % if robot_goal(2) aka y is zero. Then we must be moving straight
    v = params.max_vel;
    w = 0;
    c_input = [v; w];
    u = control_t * c_input;
elseif (R == 0) % Just rotation
    w = params.max_omega;
    v = 0;
    c_input = [v;w];
    u = control_t * c_input;
else
    if (abs(uint(1)) > abs(uint(2)))
        c = uint(1) / params.max_wheel;
    elseif (abs(uint(2)) > abs(uint(1)))
        c = uint(2) / params.max_wheel;
    elseif (abs(uint(1)) == abs(uint(2)))
        c = uint(1) / params.max_wheel;
    end
    
    u = uint / c;
end



            

        




% =================== Your code ends here ===================


end
