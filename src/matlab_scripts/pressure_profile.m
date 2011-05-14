clc;
%clear all;

P0 = 1e6;
L0 = 1.78;

% Some parameters for ease of use
U=1; V=2; W=3; RHO=5; P=4; T=6;UU=7;VV=8;WW=9;
MU=16;MUa=17;MUb=18;MUc=19;

t1 = 400;
t2 = 530;
path = '/p/lscratchd/olson45/nozzle/post_procmedium3d';
ofile = '../data/BLprof/medium.dat';

%t1 = 840;
%t2 = 880;
%path = '/p/lscratchd/olson45/nozzle/post_proccoarse3d';
%ofile = '../data/BLprof/coarse.dat'; 

if ~exist('Pmean')
    Pmean = sum_planes(path,'post',t1,t2);
    xygrid = load_planes_grid(path);
end

x = xygrid(:,:,1);
y = xygrid(:,:,2);


%plot(x(15:end,end/2)/L0,Pmean(15:end,end/2,P)/P0);
plot(x(15:end,1)/L0,Pmean(15:end,1,P)/P0);
xlim([0 8]);


%% Load the experimental mean
efile = '../data/experiment/NOZZLE_DATA/MEAN/CENTERLINE/MeanCentPress.txt'
efile = '../data/experiment/NOZZLE_DATA/MEAN/WALL/MeanWallPress.txt'
Pexp = load(efile);
figure(1);hold on;
plot(Pexp(:,1),Pexp(:,2),'ro')

