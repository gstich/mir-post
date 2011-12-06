clear all;
clc;
close all;


% Compare mean solutions to establish grid convergence.
% Toggle plots
% Figure option
LW = 2;         % LineWidth
FSn = 25;       % FontSize labels
FSa = 18;       % FontSize axis

%path = '~/Dropbox/Britton/THESIS/Figures/nozzle/convergence/';
path = '/Users/bolson/Dropbox/Britton/THESIS/Figures/nozzle/convergence/';

% Figure 1- 2d Contours of U velocity
pfig(1) = 0;
% Figure 2- Line rake of U velocity in plume
pfig(2) = 0;
% Figure 3- Line rake of U velocity in nozzle
pfig(3) = 0;
% Figure 4- Boundary layer profile, raw
pfig(4) = 0;
% Figure 5- Boundary layer profile, wall units
pfig(5) = 1;

% Load the 2 mesh sizes
load('../data/2dmv2.mat');
load('../data/2dcv2.mat');
Ht = 1.78;
Up = 32940*1.603;

% Get some variables
var_map
uc = dataC(:,:,U);
um = dataM(:,:,U);
xc = dataC(:,:,X);yc = dataC(:,:,Y);
xm = dataM(:,:,X);ym = dataM(:,:,Y);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 1- 2d Contours of U velocity %%%%%%%%
if (pfig(1) == 1);
% Contour plot
figure(1);  
cnts = linspace(-1e4,6e4,16);
cnts = linspace(-.2,1.1,16);
contour(xc/Ht,yc/Ht,uc/Up,cnts,'color','black');
hold on;
contour(xm/Ht,-ym/Ht,um/Up,cnts,'color','blue');
legend('Coarse','Medium')
xlim([-1 10]);ylim([-1 1]);
axis equal;

xlabel('X/H_t');
ylabel('Y/H_t');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 2- Line rake of U velocity in plume %
if (pfig(2) == 1);
% Line out of data in exit plume x = [0,2.5,5,7.5] Ht
Xexit = 11.7;
varC = uc; %dataC(:,:,T);
varM = um; %dataM(:,:,T);
varDIM = Up;
for i=1:4
    
xLO = Xexit + (i-1)*2.5*Ht;
% Interpolate onto line med
NN = 50;
[ZIc,XI,YIc] = lineout(xc(ceil(end/2):end,:),yc(ceil(end/2):end,:),varC(ceil(end/2):end,:),[xLO,xLO],[-3,3]*Ht,NN);
[ZIm,XI,YIm]= lineout(xm(ceil(end/2):end,:),ym(ceil(end/2):end,:),varM(ceil(end/2):end,:),[xLO,xLO],[-3,3]*Ht,NN);

figure(2);
subplot(2,2,i);
plot(ZIc/varDIM,YIc/Ht,'k-','LineWidth',3);
hold on;
plot(ZIm/varDIM,-YIm/Ht,'b-','LineWidth',3);
title(['X= ',num2str(2.5*(i-1)),'H_t'])
xlabel('U/U_p');
ylabel('Y/H_t');

ylim([-3 3]);
xlim([-.1 .6])

end
subplot(2,2,1);
legend('Coarse','Medium')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 3- Line rake of U velocity in       %
%% nozzle, through the shock
if (pfig(3) == 1);
% Line out of data in exit plume x = -[.5,1.0,1.5,2.0] Ht
Xexit = 11.7;
varC = uc; %dataC(:,:,T);
varM = um; %dataM(:,:,T);
varDIM = Up;
for i=1:4
    
xLO = Xexit - i*1*Ht;
% Interpolate onto line med
NN = 50;
[ZIc,XI,YIc] = lineout(xc(ceil(end/4):end,:),yc(ceil(end/4):end,:),varC(ceil(end/4):end,:),[xLO,xLO],[-1,1]*Ht,NN);
[ZIm,XI,YIm]= lineout(xm(ceil(end/4):end,:),ym(ceil(end/4):end,:),varM(ceil(end/4):end,:),[xLO,xLO],[-1,1]*Ht,NN);

figure(3);
subplot(2,2,i);
plot(ZIc/varDIM,YIc/Ht,'k-','LineWidth',3);
hold on;
plot(ZIm/varDIM,-YIm/Ht,'b-','LineWidth',3);
title(['X= ',num2str(-1*i),'H_t'])
xlabel('U/U_p');
ylabel('Y/H_t');

ylim([-1 1]);
xlim([-.1 1])

end
subplot(2,2,1);
legend('Coarse','Medium')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 4- BL profile at the in nozzle inlet
if (pfig(4) == 1);
[y,imC] = min( abs( xc(:,1) ));
[y,imM] = min( abs( xm(:,1) ));
    

figure(4);
yy = yc(imC,1:ceil(end/2));
yy = yy - yy(1);
uu = uc(imC,1:ceil(end/2));hold all;
%semilogx( yy / Ht, uu / Up, 'k');
plot( yy / Ht, uu / Up, 'k','LineWidth',3);

yy = ym(imM,1:ceil(end/2));
yy = yy - yy(1);
uu = um(imM,1:ceil(end/2));
%semilogx( yy / Ht, uu / Up, 'b');
plot( yy / Ht, uu / Up, 'b','LineWidth',3);
xlim([0 .15])

ylabel('U/U_p');
xlabel('Y/H_t');
legend('Coarse','Medium','Location','NorthWest')
legend boxoff;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 5- BL profile at the in nozzle inlet.. loglaw
if (pfig(5) == 1);
[y,imC] = min( abs( xc(:,1) ));
[y,imM] = min( abs( xm(:,1) ));
    

figure(5);
yy = yc(imC,1:ceil(end/2));
yy = yy - yy(1);
% Scaling
mu_w = dataC(imC,1,MU);
rho_w = dataC(imC,1,RHO);
dudy = ( dataC(imC,2,U) - dataC(imC,1,U) )/ ( yc(imC,2) - yc(imC,1) );
tauw = mu_w * dudy
utau = sqrt( tauw / rho_w );
del = mu_w / ( rho_w * utau );
uu = uc(imC,1:ceil(end/2));
semilogx( yy / del, uu / utau, 'k','LineWidth',3);hold all;

yy = ym(imM,1:ceil(end/2));
yy = yy - yy(1);
% Scaling
mu_w = dataM(imM,1,MU);
rho_w = dataM(imM,1,RHO);
dudy = ( dataM(imM,2,U) - dataM(imM,1,U) )/ ( ym(imM,2) - ym(imM,1) );
tauw = mu_w * dudy
utau = sqrt( tauw / rho_w );
del = mu_w / ( rho_w * utau );
uu = um(imM,1:ceil(end/2));
semilogx( yy / del, uu / utau, 'b','LineWidth',3);
%xlim([0 30])

x1 = logspace(-.5,1.2,20);
x2 = logspace(.8,3,20);
k = .35;
C = 5.2;
hold on;
semilogx(x1,x1,'k--');hold on;
semilogx(x2,1/k*log(x2)+C,'k--');

xlim([.5 2000])
h = ylabel('$U/U_{\tau}$'); set(h,'Interpreter','latex','FontSize',FSn);
h = xlabel('$Y/\delta_\eta$'); set(h,'Interpreter','latex','FontSize',FSn);
h = legend('Coarse','Medium','Location','NorthWest');
set(h,'Interpreter','latex','FontSize',FSn);
set(gca,'FontSize',FSa);
% Fix to make sure Latex fits in
a = get(gca,'Position');
a(1) = a(1)*1.15;
set(gca,'Position',a);

legend boxoff;

fig2pdf([path,'Ulog'],5);

end


