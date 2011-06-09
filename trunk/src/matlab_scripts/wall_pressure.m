clear all;
clc;

close all;

res = 'coarse';

% Load the entire time history of the walls
switch(res)
    case('coarse')
        path = '/p/lscratchd/olson45/nozzle/post_proccoarse3d';
        t1 = 800;%500;
        tf = 2200;%1250;
        ofile = 'DATA/Pbot_cor.mat';
    case('medium')
        path = '/p/lscratchd/olson45/nozzle/post_procmedium3d';
        t1 = 500;
        tf = 1280;%1250;
        ofile = 'DATA/Pbot_med.mat';
end




dum = load_planes_mir_xz(path,'walls_bot',t1);
Ptop = zeros(size(dum,1),size(dum,2),tf-t1+1);
Pbot = Ptop;
for i=t1:tf
  dum = load_planes_mir_xz(path,'walls_bot',i);
  Pbot(:,:,i-t1+1) = dum(:,:,1);

  %dum = load_planes_mir_xz(path,'walls_top',i);
  %Ptop(:,:,i-t1+1) = dum(:,:,1);
end

% Get xcoordinates from pressure plots
path = '/p/lscratchd/olson45/nozzle/nozzlecoarse3d/vis0000/pressure.dat';
pplot = load(path);
X = pplot(:,1);

%save 'DATA/Pbot_cor.mat' Pbot X;
save(ofile,'Pbot','X')





