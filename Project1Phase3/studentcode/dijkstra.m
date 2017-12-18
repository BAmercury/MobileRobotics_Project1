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
num_expanded = 0;

% Start and goal xy coordinates in world, convert to grid

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
%Q = zeros(length(inds)); % unvisited set, 0 means unvisited; 1 means visited

%cost = zeros(length(inds));
end_node = inds(end);

% Initialize cost
for i=1:end_node
    cost(i) = inf;
    Q(i) = 0;
    parent(i) = nan;
end

%Start node is visited
%Q(start_index) = 1;
 
% Start node cost is zero
cost(start_index) = 0;

%parent = [];

% While Q contains zeros
while ( Q(goal_index) == 0 && any( cost(:) < inf) )
    if astar == false
        lowest = inf;
        % Find node with lowest cost that is also unvisited
        for i=1:length(cost)
            b = cost(i);
            if (b < lowest)
                if (Q(i) == 0)
                    lowest = b;
                    u = i;
                end
            end
        end
    end
    %[value index] = min(cost);
    
    %disp(u);
    

    Q(u) = 1; % mark node as visited
    
    % Now we need to find neighbors of u
    [urow ucol] = ind2sub(shape,u);
    
    u_row_above = [];
    u_row_below = [];
    u_col_left = [];
    u_col_right = [];
    
    n_rows = [];
    n_cols = [];
    
    % Find neighboring rows
    if (urow == shape(1))
        u_row_above = urow - 1;
    elseif (urow == 1)
        u_row_below = urow + 1;
    else
        u_row_below = urow + 1;
        u_row_above = urow - 1;
    end
    
    % Find neighboring columns
    if (ucol == shape(2))
        u_col_left = ucol - 1;
    elseif (ucol == 1)
        u_col_right = ucol + 1;
    else
        u_col_left = ucol - 1;
        u_col_right = ucol + 1;
    end
    
    
    % Put into all in one matrix
    if ( isempty(u_row_above) == 0)
        n_rows = [n_rows, u_row_above];
    end
    n_rows = [n_rows, urow];
    if (isempty(u_row_below) == 0)
        n_rows = [n_rows, u_row_below];
    end
    

    
    if ( isempty(u_col_left) == 0)
        n_cols = [n_cols, u_col_left];
    end
    n_cols = [n_cols,ucol];
    if (isempty(u_col_right) == 0)
        n_cols = [n_cols, u_col_right];
    end
    
    neighbors = [];
    for i=1:length(n_cols)
        for j=1:length(n_rows)
            if (n_rows(j) == urow && n_cols(i) == ucol)
                neighbors = neighbors;
            else
                neighbors = [neighbors, sub2ind(shape,n_rows(j),n_cols(i))];
            end
        end
    end
%     n = [];
%     for b=1:length(neighbors)
%         v = neighbors(b);
%         if (Q(v) == 0)
%             n(i) = neighbors(b);
%         end
%     end
    
    for v=1:length(neighbors)
        
        head = grid2world(map,[urow ucol]);
        [vrow,vcol] = ind2sub(shape,neighbors(v));
        base = grid2world(map,[vrow vcol]);
        
        
        cost_between = abs( sqrt( (head(1) - (base(1)^2)) + (head(2)-(base(2)^2))   ));
        
        alt = cost(u) + cost_between;
        
        if alt < cost(neighbors(v))
            cost(neighbors(v)) = alt;
            parent(neighbors(v)) = u;
            num_expanded = num_expanded + 1;
        end
        
    end
end



path = [];
f_start = 0;
[row, col] = ind2sub(shape,goal_index);
%path = horzcat(path,[row;col]);
path = [path;row,col];
par = goal_index;
while f_start == 0
    par = parent(par);
    if par == start_index
        f_start = 1;
    end
    [row, col] = ind2sub(shape,par);
    %path = horzcat(path,[row;col]);
    path = [path;row,col];
end




% Reconstruct path




%path = [];
%for i=1:length(parent)
%    if (isnan(parent(i)) == 0)
%        [row, col] = ind2sub(shape,parent(i));
%        path = [path; row,col];
%    end
%end




    

    
path = grid2world(map,path);
    
    


    
    
%     
% path = [];
% for i=2:length(parent)
%     path = [path ;subs(parent(i),:)];
% end
% 
% path = grid2world(map,path);
% 
%     
%     
% %take p to path
% 
% %path = [];


% =================== Your code ends here ===================

end
