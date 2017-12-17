function [path, num_expanded] = dijkstra(map, start, goal, astar)
% DIJKSTRA Find the shortest path from start to goal.
%   PATH = DIJKSTRA(map, start, goal) returns an M-by-2 matrix, where each row
%   consists of the (x, y) coordinates of a point on the path.  The first
%   row is start and the last row is goal.  If no path is found, PATH is a
%   0-by-2 matrix.  Consecutive points in PATH should not be farther apart than
%   neighboring cells in the map (e.g.., if 5 consecutive points in PATH are
%   co-linear, don't simplify PATH by removing the 3 intermediate points).
%
%   PATH = DIJKSTRA(map, start, goal, astar) finds the path using euclidean
%   distance to goal as a heuristic if astar is true.
%
%   [PATH, NUM_EXPANDED] = DIJKSTRA(...) returns the path as well as
%   the number of nodes that were visited while performing the search.

% If astar is not given as an input, set the default value to be false
if nargin < 4
    astar = false;
end

% =================== Your code goes here ===================

% The example code given here may be of use to you, but it is also possible
% to solve the problem without using any of this



start = world2grid(map, start);
goal = world2grid(map, goal);


% Store the size of the map in cells
rows = map.GridSize(1);
cols = map.GridSize(2);
%n_ij = [ni, nj];
shape = [rows cols];

% Pre-compute the indices of each cell

%Returns linear index from map given [row col] aka [i j]
inds = reshape(1:rows*cols, rows, cols);

 % inds(i,j) will return the index ind of the cell with subscript [i j]
% This will return the same value as sub2ind(n_ij, i, j)

start_index = sub2ind(shape, start(1), start(2));
goal_index =  sub2ind(shape, goal(1), goal(2));

% This will return the same value as sub2ind(n_ij, i, j)

% Pre-compute the subscripts of each cell
[js, is] = meshgrid(1:cols, 1:rows);

%Given an index, get [i j], so this is the opposite of inds
subs = [is(:) js(:)];
% subs(ind,:) will return the subscript [i j] of the cell with index ind
% Note: this will return the same value as ind2sub(n_ij, ind)
% Note: subs(inds(i,j),:) will return [i j]




% Set dummy values of the outputs so that the code will run

%Algo starts here
Q = 0; % unvisited set

end_node = inds(end);

for i=1:end_node
    cost(i) = inf;
    
    %parent(i) = [];
    
    Q = [Q i];
end

%Delte that first zero
Q(1) = [];

% Start node cost is zero
cost(1) = 0;


%goal_in_Q = 1;
while ( ismember(goal_index, Q) )
    lowest = inf;
    % Find node with lowest cost within the costmap
    for i=1:length(cost)
        b = cost(i);
        if (b < lowest)
            lowest = b;
            u = i;
        end
    end
    
    Q(u) = []; % remove node from unvisisted
    
    for v=(u+1):length(Q)
        %d = cost(i) + cost(u,i) %Crunch total cost to
        space = subs(Q(v),:);
        is_occupied = getOccupancy(map,[space(1) space(2)],'grid')
        space_world = grid2world(map,[space(1) space(2)]);
        
        base = subs(Q(u),:);
        base_world = grid2world(map,[base(1) base(2)]);
        
        %if (is_occupied == 0)
            total_cost = cost(u) + abs(sqrt( ((space_world(1) - base_world(1))^2) + ((space_world(2)-space_world(2)^2))   ));
        %else
        %    total_cost = cost(u) + inf;
        %end
        
        if total_cost < cost(v)
            cost(v) = total_cost;
            parent(v) = v;
        end
        %total_cost = cost(i) + [
    end

end

path = [];
for i=2:length(parent)
    path = [path ;subs(parent(i),:)];
end

path = grid2world(map,path);

    
    
%take p to path

%path = [];
num_expanded = 0;

% =================== Your code ends here ===================

end
