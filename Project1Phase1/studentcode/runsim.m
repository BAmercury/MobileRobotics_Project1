% NOTE: This script will not run as expected unless you fill in proper
% code in trajhandle and controlhandle
% You should not modify any part of this script except for the
% visualization part
%
%% ***************** MEE 4411 SIMULATION ***************************
close all
clear all
addpath('trajectories')
addpath('utils')

% You need to implement trajhandle and controlhandle

% trajectory generator
trajhandle = @custom;

% controller
controlhandle = @controller_one;

%% *********** YOU SHOULDN'T NEED TO CHANGE ANYTHING BELOW **********

% max time
time_tol = 60;

% parameters for simulation
params = turtlebot3();

%% **************************** FIGURES *****************************
fprintf('Initializing figures...\n')
h_fig = figure;
h_3d = gca;
axis equal
grid on
xlabel('x [m]'); ylabel('y [m]')
botcolors = lines(1);

set(gcf,'Renderer','OpenGL')

%% *********************** INITIAL CONDITIONS ***********************
fprintf('Setting initial conditions...\n')
starttime = 0;          % start of simulation in seconds
tstep     = 0.01;       % this determines the time step at which the solution is given
cstep     = 0.05;       % image capture time interval
nstep     = cstep/tstep;
time      = starttime;  % current time
max_iter  = time_tol / cstep;      % max iteration
err = []; % runtime errors

% Get start and stop position
waypts    = trajhandle();
x0        = init_state(waypts(:,1), 0);
xtraj     = zeros(max_iter*nstep, length(x0));
ttraj     = zeros(max_iter*nstep, 1);

x         = x0;        % state

% Maximum position error of the robot at goal
pos_tol  = 0.02; % m
% Maximum speed of the robot at goal
vel_tol  = 0.02; % m/s

%% ************************* RUN SIMULATION *************************
fprintf('Simulation Running....')
% Main loop
for iter = 1:max_iter

    timeint = time:tstep:time+cstep;

    tic;
    % Iterate over each quad
    % Initialize quad plot
    if iter == 1
        RP = RobotPlot(x0, params.wheel_base, params.wheel_radius, 0.105, botcolors(1,:), max_iter, h_3d);
        desired_state = find_closest_point(x0(1:2), waypts);
        RP.UpdateRobotPlot(x, desired_state, time);
        h_title = title(sprintf('iteration: %d, time: %4.2f', iter, time));
    end

    % Run simulation
    [tsave, xsave] = ode45(@(t,s) botEOM(t, s, controlhandle, trajhandle, params), timeint, x);
    x              = xsave(end, :)';

    % Save to traj
    xtraj((iter-1)*nstep+1:iter*nstep,:) = xsave(1:end-1,:);
    ttraj((iter-1)*nstep+1:iter*nstep)   = tsave(1:end-1);

    % Update quad plot
    desired_state = find_closest_point(x(1:2), waypts);
    RP.UpdateRobotPlot(x, desired_state, time + cstep);

    set(h_title, 'String', sprintf('iteration: %d, time: %4.2f', iter, time + cstep))
    time = time + cstep; % Update simulation time
    t = toc;

    % Check to make sure ode45 is not timing out
    if(t > cstep*50)
        err = 'Ode45 Unstable';
        break;
    end
    % Pause to make real-time
    if (t < cstep)
        pause(cstep - t);
    end

    % Check termination criteria
    if terminate_check(x, time, waypts(:,end), pos_tol, time_tol)
        break
    end
end

%% ************************* POST PROCESSING *************************
% Truncate xtraj and ttraj
xtraj = xtraj(1:iter*nstep,:);
ttraj = ttraj(1:iter*nstep);

% Plot the saved position and velocity
% Truncate saved variables
RP.TruncateHist();
% Plot position for each quad
h_pos = figure('Name', 'Robot : position');
plot_state(h_pos, RP.state_hist(1:2,:), RP.time_hist, 'pos', 'vic');
plot_state(h_pos, RP.state_des_hist(1:2,:), RP.time_hist, 'pos', 'des');

if(~isempty(err))
    error(err);
end
fprintf('finished.\n')
