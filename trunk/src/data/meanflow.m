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

% Figure name
figs(1).name = 'u_mean_cor';
figs(2).name = 'corr_span_p';

% Plot symbols
res(1).sym = '--';
res(2).sym = '-.';
res(3).sym = '-';

%   1 U-vel 1  -4.881238E+03   4.168567E+04
%   2V-vel 1  -1.148853E+04   1.558692E+04
%   3W-vel 1  -3.861968E+03   4.124843E+03
%   4pressure 1   5.215437E+05   1.098364E+06
%   5density 1   7.277217E-04   1.462052E-03
%   6temperature 1   2.022830E+02   4.736444E+02
%   7uu 1   0.000000E+00   1.055002E+08
%   8vv 1   0.000000E+00   3.295886E+07
%   9ww 1   0.000000E+00   3.729177E+07
%   10uv 1  -3.837980E+07   2.020928E+07
%   11uw 1  -1.397939E+07   1.539493E+07
%   12vw 1  -1.035138E+07   1.123789E+07
%   13Ptot 1   5.658663E+05   1.597028E+06
%   14PtPt 1   0.000000E+00   3.453763E+10
%   gradRHO 1  -5.410950E-05   2.364354E-02
%   mu 1   6.183936E-04   1.186228E-03
%   muA 1   6.294615E-04   5.261340E-02
%   muB 1   6.268612E-04   2.042488E-02
%   muC 1   6.282928E-04   6.297713E-03

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
res(i).u = tmp(xoff:end,:,key.u);
res(i).v = tmp(xoff:end,:,key.v);
res(i).w = tmp(xoff:end,:,key.w);
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
res(i).u = tmp(xoff:end,:,key.u);
res(i).v = tmp(xoff:end,:,key.v);
res(i).w = tmp(xoff:end,:,key.w);
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


%% U-mean
vmin = -.32e4; vmax = 4e4;
Ncon = linspace(vmin,vmax,ncon);
ix = [70,250] - xoff;
iy = [32,64];
i=i+1;
figure(i);
j=1;
[cc,hh] = contourf(res(j).x,res(j).y,res(j).u,Ncon,'EdgeColor','none');
axis equal; xlim([xmin xmax]); ylim([ymin ymax]); hold on; 
cb = colorbar; h2 = ylabel(cb,'$cm/s$');set(h2,'Interpreter','latex','FontSize',FSn);
h = xlabel('$x/H_t$'); set(h,'Interpreter','latex','FontSize',FSn);
h1 = ylabel('$y/H_t$'); set(h1,'Interpreter','latex','FontSize',FSn);
plot(res(j).x(ix(1)),res(j).y(iy(1)),'ko','LineWidth',2,'MarkerFaceColor','k');hold on
plot(res(j).x(ix(2)),res(j).y(iy(2)),'ko','LineWidth',2,'MarkerFaceColor','k');
set(gca,'FontSize',FSa);


g = figure(i);
set(g,'Position',[200,395,250*AR+50,250+50])
set(g,'PaperSize',[8,8/AR]);
set(g,'PaperPositionMode','auto');

ix = [110,350] - xoff;
iy = [32,64];
i=i+1;
figure(i);
j=2;
contourf(res(j).x,res(j).y,max(res(j).u,vmin),Ncon,'EdgeColor','none');
axis equal; xlim([xmin xmax]); ylim([ymin ymax]); hold on; colorbar;
h = xlabel('$x/H_t$'); set(h,'Interpreter','latex','FontSize',FSn);
h1 = ylabel('$y/H_t$'); set(h1,'Interpreter','latex','FontSize',FSn);
plot(res(j).x(ix(1)),res(j).y(iy(1)),'ko','LineWidth',2,'MarkerFaceColor','k');hold on
plot(res(j).x(ix(2)),res(j).y(iy(2)),'ko','LineWidth',2,'MarkerFaceColor','k');




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