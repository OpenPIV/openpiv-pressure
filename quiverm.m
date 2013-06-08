function [] = quiverm(x,varargin)
% QUIVERM - plots quiver plot of matrix,
% assuming first column as X, second as Y
% third as U, and forth as V:

if isstr(x)
	x = eval(x);
end
quiver(x(:,1),x(:,2),x(:,3),x(:,4),varargin{:});

