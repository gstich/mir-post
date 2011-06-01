function data = load_planes_mir(path,field,t)

% This script will open up the 2d span averaged data and operate on themes

% Get the number of variables from the .mir file
mir_path = [path,'/mean.mir'];
%mir_path = [path,'/',field,'.mir'];
[h,d] = hdrload(mir_path);
nvar = sscanf( h(10,12:end),'%i');
len = sscanf( h(6,12:end),'%i')
nx = len(1);
ny = len(2);


% Locate the plane data and load it in raw form
data_path = [path, '/',field,num2str(t,'%04i'),'/p000000'];
disp(data_path);
funit = fopen(data_path);

   
% Put in (nx,ny) array
data = zeros(nx,ny,nvar);   
for i=1:nvar
    data(:,:,i) = fread(funit,[nx ny],'single');
end
    
   
fclose(funit);
