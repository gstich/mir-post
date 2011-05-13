clc;
%clear all;


% Some parameters for ease of use
U=1; V=2; W=3; RHO=5; P=4; T=6;UU=7;VV=8;WW=9;
MU=16;MUa=17;MUb=18;MUc=19;

t1 = 550;
t2 = 580;
path = '/p/lscratchd/olson45/nozzle/post_procmedium3d';
ofile = '../data/BLprof/medium.dat';

t1 = 800;
t2 = 840;
path = '/p/lscratchd/olson45/nozzle/post_proccoarse3d';
ofile = '../data/BLprof/coarse.dat'; 

if ~exist('Pmean')
    Pmean = sum_planes(path,'post',t1,t2);
    xygrid = load_planes_grid(path);
end

x = xygrid(:,:,1);
y = xygrid(:,:,2);


plot(x(15:end,end/2),Pmean(15:end,end/2,P));
xlim([0 20]);