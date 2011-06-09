clear all;
clc;

close all;


% Load the entire time history of the walls
path = '/p/lscratchd/olson45/nozzle/post_proccoarse3d';
t1 = 500;
tf = 1250;


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

save 'DATA/Pbot.mat' Pbot X;





