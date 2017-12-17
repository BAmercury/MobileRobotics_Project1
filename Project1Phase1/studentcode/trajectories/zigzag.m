function [waypts] = zigzag
% ZIGZAG generates waypoints to move the robot in a zigzag pattern
%
% The zigzag pattern should alternate moving the robot is the positive x
% and y directions, starting with x.  The trajectory should have 10 steps
% (a step consists of moving in x then y) and the step size (motion in x or
% y) should be 0.2m.  The robot should start at the origin.

waypts = [0; 0];
X = 0.2;
Y = 0;

for i=1:20
    waypts = [waypts, [X;Y]];
    if ( (mod(i,2) == 0)) % Check if index is even. If even, we're moving UP
        Y = Y + 0.2;
        X = X;
    else
        if (i == 1)
            Y = Y;
            X = X;
        else
            Y = Y;
            X = X+0.2;
        end
    end
    

end


end
