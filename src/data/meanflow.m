clear all;
clc;
close all


pdfE = false;

res(1).name = 'coarse';
res(2).name = 'medium';


% Figure option
LW = 2;         % LineWidth
FSn = 25;       % FontSize labels
FSa = 18;       % FontSize axis
cmax = .2;      % Contour Level Max
ncon = 32;      % Number of contour level
Ht = 1.78;
Up = 32940*1.603;
P0 = 1e6;

% Figure name
figs(1).name = 'u_mean_cor';
figs(2).name = 'u_mean_med';
figs(3).name = 'v_mean_cor';
figs(4).name = 'v_mean_med';
figs(5).name = 'p_mean_cor';
figs(6).name = 'p_mean_med';
figs(7).name = 'uu_cor';
figs(8).name = 'vv_cor';
figs(9).name = 'ww_cor';
figs(10).name = 'uv_cor';

% Plot symbols
res(1).sym = '--';
res(2).sym = '-.';
res(3).sym = '-';


key.u = 1;      key.v = 2;  key.w = 3;  key.p = 4; key.rho = 5; key.T = 6;
key.uu = 7;     key.vv = 8; key.ww = 9; key.uv = 10;
key.ptpt=14;
key.x = 20;     key.y = 21;
xoff = 5;

% Load the data
load('../matlab_scripts/DATA/MEAN_cor.mat');
load('../matlab_scripts/DATA/MEAN_med2.mat');

% Get the needed variables
vars = [key.u,key.v,key.uu,key.vv,key.ww,key.uv,key.ptpt];

% for i = max(size(vars))
%     figure(i);
%     
%     for rz = max(size(res))
%         
%         subplot( max(size(res)) , 1 , rz);

% Coarse
i=1;
tmp = MEAN_cor;
res(i).x = tmp(xoff:end,:,key.x)/Ht;
res(i).y = tmp(xoff:end,:,key.y)/Ht;
res(i).u = tmp(xoff:end,:,key.u)/Up;
res(i).v = tmp(xoff:end,:,key.v)/Up;
res(i).w = tmp(xoff:end,:,key.w)/Up;
res(i).p = tmp(xoff:end,:,key.p)/P0;
res(i).uu = tmp(xoff:end,:,key.uu);
res(i).vv = tmp(xoff:end,:,key.vv);
res(i).ww = tmp(xoff:end,:,key.ww);
res(i).uv = tmp(xoff:end,:,key.uv);
res(i).ptpt = tmp(xoff:end,:,key.ptpt);
clear tmp MEAN_cor;

% Medium
i=2;
tmp = MEAN_med;
res(i).x = tmp(xoff:end,:,key.x)/Ht;
res(i).y = tmp(xoff:end,:,key.y)/Ht;
res(i).u = tmp(xoff:end,:,key.u)/Up;
res(i).v = tmp(xoff:end,:,key.v)/Up;
res(i).w = tmp(xoff:end,:,key.w)/Up;
res(i).p = tmp(xoff:end,:,key.p)/P0;
res(i).uu = tmp(xoff:end,:,key.uu);
res(i).vv = tmp(xoff:end,:,key.vv);
res(i).ww = tmp(xoff:end,:,key.ww);
res(i).uv = tmp(xoff:end,:,key.uv);
res(i).ptpt = tmp(xoff:end,:,key.ptpt);
clear tmp MEAN_med;


xmin = 0;
xmax = 8;
ymin = -1;
ymax = 1;

AR = (xmax-xmin) / (ymax-ymin);

i=0;
xrng = [xmin xmax];
yrng = [ymin ymax];

j=1;iFig=1;
vrng = [-.1,.8];
xL = '$x/H_t$';
yL = '$y/H_t$';
cbL = '$u/U_p$';
mean_contour(res(j).x,res(j).y,res(j).u,xrng,yrng,vrng,ncon,iFig,xL,yL,cbL);
ix = [380,87] - xoff;
iy = [80,18];
plot(res(j).x(ix(1),iy(1)),res(j).y(ix(1),iy(1)),'go','LineWidth',2,'MarkerFaceColor','k','MarkerSize',10);hold on
plot(res(j).x(ix(2),iy(2)),res(j).y(ix(2),iy(2)),'go','LineWidth',2,'MarkerFaceColor','k','MarkerSize',10);

j=2;iFig=iFig+1;
vrng = [-.1,.8];
xL = '$x/H_t$';
yL = '$y/H_t$';
cbL = '$u/U_p$';
mean_contour(res(j).x,res(j).y,res(j).u,xrng,yrng,vrng,ncon,iFig,xL,yL,cbL);
ix = [577,131] - xoff;
iy = [95,21];
plot(res(j).x(ix(1),iy(1)),res(j).y(ix(1),iy(1)),'go','LineWidth',2,'MarkerFaceColor','k','MarkerSize',10);hold on
plot(res(j).x(ix(2),iy(2)),res(j).y(ix(2),iy(2)),'go','LineWidth',2,'MarkerFaceColor','k','MarkerSize',10);


%% V-velocity
j=1;iFig=iFig+1;
vrng = [-.1,.1];
xL = '$x/H_t$';
yL = '$y/H_t$';
cbL = '$v/U_p$';
mean_contour(res(j).x,res(j).y,res(j).v,xrng,yrng,vrng,ncon,iFig,xL,yL,cbL);

j=2;iFig=iFig+1;
vrng = [-.1,.1];
xL = '$x/H_t$';
yL = '$y/H_t$';
cbL = '$v/U_p$';
mean_contour(res(j).x,res(j).y,res(j).v,xrng,yrng,vrng,ncon,iFig,xL,yL,cbL);


%% Pressure
j=1;iFig=iFig+1;
vrng = [.4,1];
xL = '$x/H_t$';
yL = '$y/H_t$';
cbL = '$p/p_0$';
mean_contour(res(j).x,res(j).y,res(j).p,xrng,yrng,vrng,ncon,iFig,xL,yL,cbL);

j=2;iFig=iFig+1;
vrng = [.4,1];
xL = '$x/H_t$';
yL = '$y/H_t$';
cbL = '$p/p_0$';
mean_contour(res(j).x,res(j).y,res(j).p,xrng,yrng,vrng,ncon,iFig,xL,yL,cbL);

%% u'u'
j=1;iFig=iFig+1;
vrng = [0,0];
xL = '$x/H_t$';
yL = '$y/H_t$';
cbL = '$u_{rms}$';
mean_contour(res(j).x,res(j).y,res(j).uu,xrng,yrng,vrng,ncon,iFig,xL,yL,cbL);

%% v'v'
j=1;iFig=iFig+1;
vrng = [0,0];
xL = '$x/H_t$';
yL = '$y/H_t$';
cbL = '$v_{rms}$';
mean_contour(res(j).x,res(j).y,res(j).vv,xrng,yrng,vrng,ncon,iFig,xL,yL,cbL);

%% w'w'
j=1;iFig=iFig+1;
vrng = [0,0];
xL = '$x/H_t$';
yL = '$y/H_t$';
cbL = '$w_{rms}$';
mean_contour(res(j).x,res(j).y,res(j).ww,xrng,yrng,vrng,ncon,iFig,xL,yL,cbL);

%% u'v''
j=1;iFig=iFig+1;
vrng = [0,0];
xL = '$x/H_t$';
yL = '$y/H_t$';
cbL = '$uv_{rms}$';
mean_contour(res(j).x,res(j).y,res(j).uv,xrng,yrng,vrng,ncon,iFig,xL,yL,cbL);


%% Ptot'Ptot'
j=1;iFig=iFig+1;
vrng = [0,0];
xL = '$x/H_t$';
yL = '$y/H_t$';
cbL = '$Pt_{rms}$';
mean_contour(res(j).x,res(j).y,res(j).ptpt,xrng,yrng,vrng,64,iFig,xL,yL,cbL);




% Save the figures and convert them to .pdf
if (pdfE)
    for i=1 : size (figs , 2)
        fname = [ '../figs/',figs(i).name , '.eps' ];
        figure(i);
        print('-depsc2',fname)
        eps2pdf(fname)
        delete(fname)
    end
end