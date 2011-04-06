function data = load_planes(path,field,t)

% This script will open up the 2d span averaged data and operate on themes

% Locate the plane data and load it in raw form
data_path = [path, '/vis',num2str(t,'%04i'),'/planes/',field,'.tec'];
disp(data_path);
[h,d] = hdrload(data_path);

% Parse th nx and ny
eqs = findstr(h(2,:),'=');
comma = findstr(h(2,:),',');
nx = str2double( h(2,(eqs(1)+1:comma(1)-1)));
ny = str2double( h(2,(eqs(2)+1:comma(2)-1)));
nvar = 8;
   
% Put in (nx,ny) array
data = zeros(nx,ny,nvar);   
count = 1;
for j=1:ny
    for i=1:nx
        data(i,j,1:2) = d(count  ,1:2);
        data(i,j,3:4) = d(count+1,1:2);
        data(i,j,5:6) = d(count+2,1:2);
        data(i,j,7:8) = d(count+3,1:2);
        count = count + 4;
    end
end
    
   
