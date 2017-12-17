function [x] = qdToState(qd)
% Converts state vector for simulation to qd struct used in hardware.
% x is 1 x 3 vector of state variables [pos theta]
% qd is a struct including the fields pos, theta

x = zeros(1,3); %initialize dimensions

x(1:2) = qd.pos;
x(3) = qd.theta;

end
