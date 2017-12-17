function [theta] = RotToTheta(R)
%RotToTheta Extract angle from a world-to-body Rotation Matrix
%   The rotation matrix in this function is world to body [bRw] you will
%   need to transpose the matrix if you have a body to world [wRb] such
%   that [wP] = [wRb] * [bP], where [bP] is a point in the body frame and
%   [wP] is a point in the world frame

phi = atan2(R(2,1), R(1,1));

end
