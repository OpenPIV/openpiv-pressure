function [x,y,u,v] = parse_tmp_matrix(tmp,size_x,size_y)
% [x,y,u,v] = parse_tmp_matrix(tmp,size_x,size_y)
% assigns the columns in the 'tmp' matrix according to the order 
% of the parameters: 1-x,2-y,3-u,4-v and reshapes them in rectangular
% matrices according to the provided input of 'size_x','size_y'
% no error checks at the moment 

x = reshape(tmp(:,1),size_x,size_y);
y = reshape(tmp(:,2),size_x,size_y);
u = reshape(tmp(:,3),size_x,size_y);
v = reshape(tmp(:,4),size_x,size_y);