function R = ThetaToRot(theta)
%ThetaToRot Converts angle to a body-to-world Rotation matrix
%   The rotation matrix in this function is world to body [bRw] you will
%   need to transpose this matrix to get the body to world [wRb] such that
%   [wP] = [wRb] * [bP], where [bP] is a point in the body frame and [wP]
%   is a point in the world frame
%   written by Daniel Mellinger
%

R = [cos(theta), -sin(theta);
     sin(theta), cos(theta)];

end
