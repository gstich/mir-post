function data = sum_planes(path,field,t1,t2)

% This script will averaged a series of planes and return the resultant
% plane
data = 0;
for tt=t1:t2
    plane = load_planes(path,field,tt);
    data = data + plane/ (t2-t1+1);
end
   


