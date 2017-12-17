function [qd] = stateToQd(x)
%Converts qd struct used in hardware to x vector used in simulation
% x is 1 x 3 vector of state variables [pos theta]
% qd is a struct including the fields pos, theta

%current state
qd.pos = x(1:2);
qd.theta = x(3);

end
