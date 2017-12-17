function map = denvToOccGrid(boundary, blocks, sz, margin)
% ENVTOOCCGRID
%
% This function takes in the output of readEnv and two additional
% parameters and outputs BinaryOccupancyGrid object that represents the
% input environment.
%
% INPUTS:
% boundary  - 1 x 4, matrix with [xmin, ymin, xmax, ymax] of the environment
% blocks 	- n x 4, matrix with n rows, one for each block, of the form
%                    [xmin, ymin, xmax, ymax]
% sz        - 1 x 1, size of the cells [m/cell]
% margin    - 1 x 1, amount to inflate each obstacle by (in both x and y)
%
% OUTPUTS:
% map       - robotics.BinaryOccupancyGrid object

% =================== Your code goes here ===================

% This opens the documentation for the occupancy grid. You should comment 
% this out later if you don't want it popping up every time you run the 
% function.
%doc robotics.BinaryOccupancyGrid 

% You should fill this in!  I picked dummy values just so it would run.
width = norm(boundary(1) - boundary(3));
height = norm(boundary(2) - boundary(4));
map = robotics.BinaryOccupancyGrid(width,height,(1/sz));

map.GridLocationInWorld = [boundary(1) boundary(2)]

%% Find out what cells overlap grid

% Separate blocks into a matrix of: [xmin ymin; xmax ymax] instead of [xmin
% ymin xmax ymax; xmin ymin xmax ymax]
block_size = size(blocks);
block_rows = block_size(1);
points = [];
for g=1:block_rows
    points = [points; blocks(g,1) blocks(g,2);blocks(g,3) blocks(g,4)];
end

%points = [blocks(:,1) blocks(:,2); blocks(:,3) blocks(:,4)];
%points = [blocks(1) blocks(4); blocks(2) blocks(5); blocks(3) blocks(6); blocks(7) blocks(10); blocks(8) blocks(11); blocks(9) blocks(12)]

grid_points = world2grid(map, points);



for k=1:2:(length(grid_points)-1)

    grid_point_block = [grid_points(k,:); grid_points(k+1,:)];
    
    % Find min_row, max_row, min_col, and max_col
    %rows:
    
    if (grid_point_block(1,1) > grid_point_block(2,1))
        max_row = grid_point_block(1,1);
        min_row = grid_point_block(2,1);
    else
        max_row = grid_point_block(2,1);
        min_row = grid_point_block(1,1);
    end
    
        if (grid_point_block(1,2) > grid_point_block(2,2))
        max_col = grid_point_block(1,2);
        min_col = grid_point_block(2,2);
    else
        max_col = grid_point_block(2,2);
        min_col = grid_point_block(1,2);
    end


    
    for (i=min_row:max_row)
        row_val = ismember(i,grid_point_block(:,1));

        for (j=min_col:max_col)
           col_val = ismember(j,grid_point_block(:,2));

           setOccupancy(map,[i j],1,'grid')

        end


end
    
    
    




map.inflate(margin)
% Width, Height, Res

% =================== Your code ends here ===================
end
