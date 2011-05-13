function show_points_cor
clc;
close all;

% Some parameters for ease of use
U=1; V=2; W=3; RHO=5; P=4; T=6;UU=7;VV=8;WW=9;
MU=16;MUa=17;MUb=18;MUc=19;

VAR = U;
Ncon = 32;
ny = 128;
yoff = 25;

t1 = 800;
t2 = 840;
x1 = 200;
xn = 440;
path = '/p/lscratchd/olson45/nozzle/post_proccoarse3d';
ifile = '~/miranda/chaos_4_x86_64_ib/bin/corr.cor';

if ~exist('Pmean')
    Pmean = sum_planes(path,'post',t1,t2);
end

grid2d = load_planes_grid(path);
X = grid2d(x1:xn,:,1);
Y = grid2d(x1:xn,:,2);
V = Pmean(x1:xn,:,VAR);


contourf(X,Y,V,Ncon,'edgecolor','none');
%plot(

% Load the file with correlation points
%Pcorr = load(ifile);

Pcorr = grid_noz(240,ny-2*yoff,8*3,7*3);
for i=2:Pcorr(1,1)+1
    ix = Pcorr(i,1);
    iy = Pcorr(i,2)+yoff;
    hold on; plot(X(ix,iy),Y(ix,iy),'ko')
end

Pcorr(2:end,1) = Pcorr(2:end,1) + x1;
Pcorr(2:end,2) = Pcorr(2:end,2) + yoff;
dlmwrite(ifile,Pcorr,'delimiter',' ');

end



function Pcorr = grid_noz(nx,ny,mx,my)


ix = linspace(1,nx,mx);
ix = ceil(ix);

iy = linspace(1,ny,my);
iy = ceil(iy);

Pcorr(1:mx*my,1:2)=0;
Pcorr(1,1) = mx*my;
Pcorr(1,2) = mx;

for i=2:mx*my+1
    
end

ii=1;
for j=1:my
    for i=1:mx
        ii = ii + 1;
        Pcorr(ii,1) = ix(i);
        Pcorr(ii,2) = iy(j);
    end
end






end

