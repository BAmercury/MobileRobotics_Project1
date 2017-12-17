classdef RobotPlot < handle
    %BOTPLOT Visualization class for ground robot

    properties (SetAccess = public)
        k = 0;
        qn;             % robot number
        time = 0;       % time
        state;          % state
        des_state;      % desired state [x; y];
        rot;            % rotation matrix body to world

        color;          % color of robot
        wheel_base;     % wheelbase
        height;         % height of robt
        wheel_radius;   % wheel radius
        radius;         % radius of robot
        wheel;          % wheel positions

        state_hist      % position history
        state_des_hist; % desired position history
        time_hist;      % time history
        max_iter;       % max iteration
    end

    properties (SetAccess = private)
        h_2d
        h_lw;           % left wheel handle
        h_rw;           % right wheel handle
        h_base;         % robot base handle
        h_qn;           % robot number handle
        h_pos_hist;     % position history handle
        h_pos_des_hist; % desired position history handle
        text_dist;      % distance of robot number to robot
    end

    methods
        % Constructor
        function Q = RobotPlot(state, wheel_base, wheel_radius, radius, color, max_iter, h_2d)
            Q.state = state(1:2);
            Q.wheel_base = wheel_base;
            Q.color = color;
            Q.wheel_radius = wheel_radius;
            Q.radius = radius;
            Q.rot = ThetaToRot(state(3));
            Q.text_dist = Q.radius / 3;
            Q.des_state = Q.state(1:2);
            Q.wheel = robot_pos(Q.state(1:2), Q.rot, Q.wheel_radius, Q.wheel_base, Q.radius);

            Q.max_iter = max_iter;
            Q.state_hist = zeros(2, max_iter);
            Q.state_des_hist = zeros(2, max_iter);
            Q.time_hist = zeros(1, max_iter);

            % Initialize plot handle
            if nargin < 7, h_2d = gca; end
            Q.h_2d = h_2d;
            hold(Q.h_2d, 'on')
            Q.h_pos_hist = plot(Q.h_2d, Q.state(1), Q.state(2), 'r.');
            Q.h_pos_des_hist = plot(Q.h_2d, Q.des_state(1), Q.des_state(2), 'b.');
            Q.h_base = fill(Q.h_2d, ...
                Q.wheel(1,5:end), ...
                Q.wheel(2,5:end), ...
                'k', ...
                'FaceColor', Q.color, 'LineWidth', 2);
            Q.h_lw = plot(Q.h_2d, ...
                Q.wheel(1,[1 2]), ...
                Q.wheel(2,[1 2]), ...
                '-k', 'MarkerFaceColor', Q.color, 'MarkerSize', 5);
            Q.h_rw = plot(Q.h_2d, ...
                Q.wheel(1,[3 4]), ...
                Q.wheel(2,[3 4]), ...
                '-k', 'MarkerFaceColor', Q.color, 'MarkerSize', 5);
            hold(Q.h_2d, 'off')
        end

        % Update quad state
        function UpdateRobotState(Q, state, time)
            Q.state = state;
            Q.time = time;
            Q.rot = ThetaToRot(state(3)); % Q.rot needs to be body-to-world
        end

        % Update desired quad state
        function UpdateDesiredRobotState(Q, des_state)
            Q.des_state = des_state;
        end

        % Update quad history
        function UpdateRobotHist(Q)
            Q.k = Q.k + 1;
            Q.time_hist(Q.k) = Q.time;
            Q.state_hist(:,Q.k) = Q.state(1:2);
            Q.state_des_hist(:,Q.k) = Q.des_state(1:2);
        end

        % Update motor position
        function UpdateWheelPos(Q)
            Q.wheel = robot_pos(Q.state(1:2), Q.rot, Q.wheel_radius, Q.wheel_base, Q.radius);
        end

        % Truncate history
        function TruncateHist(Q)
            Q.time_hist = Q.time_hist(1:Q.k);
            Q.state_hist = Q.state_hist(:, 1:Q.k);
            Q.state_des_hist = Q.state_des_hist(:, 1:Q.k);
        end

        % Update quad plot
        function UpdateRobotPlot(Q, state, des_state, time)
            Q.UpdateRobotState(state, time);
            Q.UpdateDesiredRobotState(des_state);
            Q.UpdateRobotHist();
            Q.UpdateWheelPos();
            set(Q.h_lw, ...
                'XData', Q.wheel(1,[1 2]), ...
                'YData', Q.wheel(2,[1 2]));
            set(Q.h_rw, ...
                'XData', Q.wheel(1,[3 4]), ...
                'YData', Q.wheel(2,[3 4]));
            set(Q.h_base, ...
                'XData', Q.wheel(1,5:end), ...
                'YData', Q.wheel(2,5:end))
            set(Q.h_pos_hist, ...
                'XData', Q.state_hist(1,1:Q.k), ...
                'YData', Q.state_hist(2,1:Q.k));
            set(Q.h_pos_des_hist, ...
                'XData', Q.state_des_hist(1,1:Q.k), ...
                'YData', Q.state_des_hist(2,1:Q.k));
            drawnow;
        end
    end

end
