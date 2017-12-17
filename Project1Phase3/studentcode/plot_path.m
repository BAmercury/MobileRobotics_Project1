function plot_path(map, path)
% PLOT_PATH Visualize a path through an environment
%   PLOT_PATH(map, path) creates a figure showing a path through the
%   environment.  path is an N-by-2 matrix where each row corresponds to the
%   (x, y) coordinates of one point along the path.

% Plot the map
map.show;

% Plot the path
hold on;
plot(path(:,1), path(:,2), 'r-');
hold off;

end