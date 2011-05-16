clc;
%clear all;

P0 = 1e6;
L0 = 1.78;

% Some parameters for ease of use
U=1; V=2; W=3; RHO=5; P=4; T=6;UU=7;VV=8;WW=9;
MU=16;MUa=17;MUb=18;MUc=19;

resolution = 'medium';
%resolution = 'coarse';

switch resolution

    case 'medium'
    t1 = 400;
    t2 = 530;
    path = '/p/lscratchd/olson45/nozzle/post_procmedium3d';
    ofile = '../data/pressure/medium.dat';
    xoff = -.24;

    case 'coarse'
    t1 = 1000;
    t2 = 1200;
    path = '/p/lscratchd/olson45/nozzle/post_proccoarse3d';
    ofile = '../data/pressure/coarse.dat'; 
    xoff = -.24;

end
if ~exist('Pmean')
    Pmean = sum_planes(path,'post',t1,t2);
    xygrid = load_planes_grid(path);
end

x = xygrid(:,:,1);
y = xygrid(:,:,2);

% Centerline
figure(1);
plot(x(15:end,end/2)/L0 + xoff,Pmean(15:end,end/2,P)/P0);
xlim([0 8]);
% Wall
figure(2);hold on;
plot(x(15:end,end)/L0 + xoff,Pmean(15:end,end,P)/P0);
plot(x(15:end,end)/L0 + xoff,Pmean(15:end,1,P)/P0);
xlim([0 8]);


%% Load the experimental mean
efile = '../data/experiment/NOZZLE_DATA/MEAN/CENTERLINE/MeanCentPress.txt'
Pexp = load(efile);
figure(1);hold on;
plot(Pexp(:,1),Pexp(:,2),'ro')

efile = '../data/experiment/NOZZLE_DATA/MEAN/WALL/MeanWallPress.txt'
Pexp = load(efile);
figure(2);hold on;
plot(Pexp(:,1),Pexp(:,2),'ro')

clear press;
press(:,1) = x(:,end/2)/L0;
press(:,2) = Pmean(:,end/2,P)/P0;
press(:,3) = x(:,1)/L0;
press(:,4) = Pmean(:,1,P)/P0;
press(:,5) = Pmean(:,end,P)/P0;

save ofile press
key1 = '%% < 1-4  > x/H, P_cen/P0, x/H, P_wall/P0 (lower),(upper)';
dlmwrite(ofile,key1,'delimiter',' ');
dlmwrite(ofile,press,'delimiter',' ','-append');


