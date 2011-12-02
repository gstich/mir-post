function data = load_planes_grid(path)

% This script will open up the 2d span averaged data and operate on themes

% Get the number of variables from the .mir file
mir_path = [path,'/mean.mir'];
mir_path = [path,'/post.mir'];
[h,d] = hdrload(mir_path);
len = sscanf( h(6,12:end),'%i');
nx = len(1);
ny = len(2);


% Locate the plane data and load it in raw form
data_path = [path, '/grid2d/p000000'];
disp(data_path);
funit = fopen(data_path);

   
% Put in (nx,ny) array
data = zeros(nx,ny,2);   
for i=1:2
    data(:,:,i) = fread(funit,[nx ny],'single');
end
    
fclose(funit);
