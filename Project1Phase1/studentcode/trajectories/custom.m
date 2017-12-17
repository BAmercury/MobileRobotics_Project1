function [waypts] = custom



waypts = pathFind2(10,300,400);

waypts = waypts .* 0.01;

waypts = transpose(waypts);

%waypts = [0, waypts(1,:); 0, waypts(2,:)]








end
