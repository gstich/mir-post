function data = sum_planes(path,field,t1,t2)

% This script will averaged a series of planes and return the resultant
% plane
data = 0;
for tt=t1:t2
    %plane = load_planes(path,field,tt);
    plane = load_planes_mir(path,field,tt);
    data = data + plane/ (t2-t1+1);
end
   

% Get the grid and tak in on to the data array
grid = load_planes_grid(path);

nvar = size(data,3);
data(:,:,nvar+1:nvar+2) = grid;

% Fix the wierd thing at the inlet.. chop off data there.
bdata = data;
clear data;

data = bdata(10:end,:,:);
