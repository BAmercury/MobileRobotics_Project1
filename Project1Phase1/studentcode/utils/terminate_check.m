function [ terminate_cond ] = terminate_check(x, time, stop, pos_tol, time_tol)
%TERMINATE_CHECK Check termination criteria, including position, velocity and time

if time < 5
    terminate_cond = 0;
    return;
end

% Initialize
pos_check = true;

% Check position and velocity and still time for each quad
pos_check = pos_check && (norm(x(1:2) - stop) < pos_tol);

% Check total simulation time
time_check = time > time_tol;

if pos_check
    terminate_cond = 1; % Robot reaches goal and stops, successful
elseif time_check
    terminate_cond = 2; % Robot doesn't reach goal within given time, not complete
else
    terminate_cond = 0;
end

end