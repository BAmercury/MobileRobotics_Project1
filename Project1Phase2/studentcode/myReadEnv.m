function [boundary, blocks] = myReadEnv(filename)
% MYREADENV
%
% This function takes in the name of a text file and creates two matrices,
% one that describes the size of the environment and one that describes the
% size and locations of obstacles in the environment.  The map and the
% obstacles are assumed to be axis-aligned rectangles.  Example inputs are
%   boundary 1.0 2.0 4.0 5.0
% which means that the lower left corner of the environment is at (1,2) and
% the upper right is at (4,5).
%
% INPUTS:
% filename 	- string giving the name and location of the input file
%
% OUTPUTS:
% boundary  - 1 x 4, matrix with [xmin, ymin, xmax, ymax] of the environment
% blocks 	- n x 4, matrix with n rows, one for each block, of the form
%                    [xmin, ymin, xmax, ymax]

% =================== Your code goes here ===================

% You should fill this in!  Set to empty matrix so that the function will
% run when you call it.


data = importdata(filename);

boundary = data.data(1,:);

blocks = data.data(2:end,:);

% =================== Your code ends here ===================

end
