clear all;
clc;

close all;

res = 'medium';

% Load the entire time history of the walls
switch(res)
    case('coarse')
        path = '/p/lscratchd/olson45/nozzle/post_proccoarse3d';
        t1 = 800;%500;
        tf = 2258;%1250;
        ofile = 'DATA/Pw_cor.mat';
        side = 'walls_top';
    case('medium')
        path = '/p/lscratchd/olson45/nozzle/post_procmedium3d';
        t1 = 300;
        tf = 1280;%1250;
        ofile = 'DATA/Pw_med.mat';
        side = 'walls_bot';        
end




dum = load_planes_mir_xz(path,'walls_bot',t1);
Pw = zeros(size(dum,1),size(dum,2),tf-t1+1);
for i=t1:tf
  dum = load_planes_mir_xz(path,side,i);
  Pw(:,:,i-t1+1) = dum(:,:,1);

end

% Get xcoordinates from pressure plots
path = ['/p/lscratchd/olson45/nozzle/nozzle',res,'3d/vis0000/pressure.dat'];
pplot = load(path);
X = pplot(:,1);

%save 'DATA/Pbot_cor.mat' Pbot X;
save(ofile,'Pw','X')





