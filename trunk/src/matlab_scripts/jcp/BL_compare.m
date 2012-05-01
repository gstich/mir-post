clear all;
clc;
close all;


path(path,'../');
% Compare mean solutions to establish grid convergence.
% Toggle plots
% Figure option
LW = 2;         % LineWidth
FSn = 25;       % FontSize labels
FSa = 18;       % FontSize axis

pdfE = false;
%pdfE = true;

%path = '~/Dropbox/Britton/THESIS/Figures/nozzle/convergence/';
path = '/Users/bolson/Dropbox/Britton/THESIS/Figures/nozzle/convergence/';

% Figure 1- 2d Contours of U velocity
pfig(1) = 0;
f(1).name = 'Ucontours';

% Figure 4- Boundary layer profile, raw
pfig(4) = 1;
f(4).name = 'UBLmean';

% Figure 5- Boundary layer profile, wall units
pfig(5) = 1;
f(5).name = 'UBLlog';

% Figure 6- Artificial terms
pfig(6) = 0;
f(6).name = 'artMU';

% Figure 11- BL RMS profiles
pfig(11) = 1;
f(11).name = 'BL_RMS';




%pfig(:) = 0;

var_map

% Load the 2 mesh sizes
load('../../data/jcp/fin_muNew.mat');
dataM = data;

load('../../data/jcp/fin_nomu.mat');
dataF = data;

Ht = 0.2;
Up = data(1,end,U)
Pamb = 1e6;
delBL = 0.2;

% Get some variables

um = sqrt( dataM(:,:,U).^2 + dataM(:,:,V).^2 );
uf = sqrt( dataF(:,:,U).^2 + dataF(:,:,V).^2 );

xm = dataM(:,:,X);ym = dataM(:,:,Y);
xf = dataF(:,:,X);yf = dataF(:,:,Y);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 4- BL profile at the in nozzle inlet
if (pfig(4) == 1);

[y,imM] = min( abs( xm(:,1)-2.5 ));
[y,imF] = min( abs( xf(:,1)-2.5 ));
    

figure(4);
yy = ym(imM,1:ceil(end/2));
yy = yy - yy(1);
uu = um(imM,1:ceil(end/2));hold all
%semilogx( yy / Ht, uu / Up, 'b');
plot( yy / Ht, uu / Up, 'b','LineWidth',LW);

yy = yf(imF,1:ceil(end/2));
yy = yy - yy(1);
uu = uf(imF,1:ceil(end/2));
%semilogx( yy / Ht, uu / Up, 'b');
plot( yy / Ht, uu / Up, 'k--','LineWidth',LW);

xlim([0 .15])

h = ylabel('$U/U_p$'); set(h,'Interpreter','latex','FontSize',FSn);
h = xlabel('$Y/H_t$'); set(h,'Interpreter','latex','FontSize',FSn);
h = legend('Mesh A','Mesh B','Location','NorthWest');
set(h,'Interpreter','latex','FontSize',FSn);
legend boxoff;
box on;

if (pdfE)
fig2pdf([path,f(4).name],4);
end


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 5- BL profile at the in nozzle inlet.. loglaw
if (pfig(5) == 1);

[y,imM] = min( abs( xm(:,1)-2.5 ));
[y,imF] = min( abs( xf(:,1)-2.5 ));
    

figure(5);

yy = ym(imM,1:ceil(end/2));
yy = yy - yy(1);
% Scaling
mu_w = dataM(imM,1,MU);
rho_w = dataM(imM,1,RHO);
dudy = ( um(imM,2) - um(imM,1) )/ ( ym(imM,2) - ym(imM,1) );
tauw = mu_w * dudy
utau = sqrt( tauw / rho_w );
del = mu_w / ( rho_w * utau );
uu = um(imM,1:ceil(end/2));

% Van Driest Transform
uvd = zeros(1,max(size(uu)));
for j=1:max(size(uu))
  uvd(j) = trapz( uu(1:j),sqrt( dataM(imM,1:j,RHO)/rho_w ) );
end 



semilogx( yy / del, uvd / utau, 'bo','LineWidth',LW);hold all;

yy = yf(imF,1:ceil(end/2));
yy = yy - yy(1);
% Scaling
mu_w = dataF(imF,1,MU);
rho_w = dataF(imF,1,RHO);
dudy = ( uf(imF,2) - uf(imF,1) )/ ( yf(imF,2) - yf(imF,1) );
tauw = mu_w * dudy
utau = sqrt( tauw / rho_w );
del = mu_w / ( rho_w * utau );
uu = uf(imF,1:ceil(end/2));

% Van Driest
%uvd(1) = 0;
%for i=2:max(size(uu))
%    dup = uu(i) - uu(i-1);
%    uvd(i) = uvd(i-1) + sqrt( dataF(imF,i,RHO) / rho_w) * dup;
%end
% Van Driest Transform
uvd = zeros(1,max(size(uu)));
for j=1:max(size(uu))
  uvd(j) = trapz( uu(1:j),sqrt( dataF(imF,1:j,RHO)/rho_w ) );
end 


semilogx( yy / del, uvd / utau, 'k--','LineWidth',LW);


%xlim([0 30])

x1 = logspace(-.5,1.2,20);
x2 = logspace(.8,3,20);
k = .41;
C = 5.1;
hold on;
semilogx(x1,x1,'k--');hold on;
semilogx(x2,1/k*log(x2)+C,'k--');

xlim([.5 2000])
h = ylabel('$u_{vd}^+$'); set(h,'Interpreter','latex','FontSize',FSn);
h = xlabel('$y^+$'); set(h,'Interpreter','latex','FontSize',FSn);
h = legend('Mesh A','Mesh B','Location','NorthWest');
set(h,'Interpreter','latex','FontSize',FSn);
set(gca,'FontSize',FSa);
% Fix to make sure Latex fits in
a = get(gca,'Position');
a(1) = a(1)*1.15;
set(gca,'Position',a);

legend boxoff;

if (pdfE)
fig2pdf([path,f(5).name],5);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 6- 
if (pfig(6) == 1);
figure(6);hold all;
    
[vf2,xif,yif] = lineout(xf,yf,dataF(:,:,MU),[2.7,2.7],[0,0.2],50);
[vm2,xim,yim] = lineout(xm,ym,dataM(:,:,MU),[2.7,2.7],[0,0.2],50);
display('1 of 3');

[vf,xif,yif] = lineout(xf,yf,dataF(:,:,MUb),[2.7,2.7],[0,0.2],50);
[vm,xim,yim] = lineout(xm,ym,dataM(:,:,MUb),[2.7,2.7],[0,0.2],50);
display('2 of 3');

plot(yim/Ht,vm./vm2,'b','Linewidth',LW);
plot(yif/Ht,vf./vf2,'r','Linewidth',LW);

[vf,xif,yif] = lineout(xf,yf,dataF(:,:,MUa),[11.7,11.7],[-1.347,1.347],50);
[vm,xim,yim] = lineout(xm,ym,dataM(:,:,MUa),[11.7,11.7],[-1.347,1.347],50);
display('3 of 3');

plot(yim/Ht,vm./vm2,'b--','Linewidth',LW/2);
plot(yif/Ht,vf./vf2,'r--','Linewidth',LW);

h = legend('$\mu^*_{\eta}$ Mesh A','$\mu^*_{\eta}$ Mesh B','$\mu^*_{\eta}$ Mesh C');
set(h,'Interpreter','latex','FontSize',FSn);
legend boxoff;
box on;

h = xlabel('$Y/H_t$'); set(h,'Interpreter','latex','FontSize',FSn);
h = ylabel('$\mu*/\mu$'); set(h,'Interpreter','latex','FontSize',FSn);

if (pdfE)
fig2pdf([path,f(6).name],6);
end


end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 11- BL profile at the in nozzle inlet
if (pfig(11) == 1);

[y,imM] = min( abs( xm(:,1)-2.5 ));
[y,imF] = min( abs( xf(:,1)-2.5 ));
    
ivar = [UU,VV,WW];

figure(11);


% Scaling
mu_w = dataM(imM,1,MU);
rho_w = dataM(imM,1,RHO);
dudy = ( um(imM,2) - um(imM,1) )/ ( ym(imM,2) - ym(imM,1) );
tauw = mu_w * dudy
utau = sqrt( tauw / rho_w );
del = mu_w / ( rho_w * utau );
cf = 2*tauw / dataM(imM,end,RHO)/um(imM,end)^2

for iv = 1:max(size(ivar));
yy = ym(imM,1:ceil(end/2));
yy = yy - yy(1);
uu = dataM(imM,1:ceil(end/2),ivar(iv));
uu = sqrt(uu) ./ utau .*sqrt(dataM(imM,1:ceil(end/2),RHO)/rho_w);
hold all;
plot( yy / delBL, uu, 'b','LineWidth',LW);
end

% Scaling
mu_w = dataF(imF,1,MU);
rho_w = dataF(imF,1,RHO);
dudy = ( uf(imF,2) - uf(imF,1) )/ ( yf(imF,2) - yf(imF,1) );
tauw = mu_w * dudy
utau = sqrt( tauw / rho_w );
del = mu_w / ( rho_w * utau );
cf = 2*tauw / dataF(imF,end,RHO)/uf(imF,end)^2

for iv = 1:max(size(ivar));
yy = yf(imF,1:ceil(end/2));
yy = yy - yy(1);
uu = dataF(imF,1:ceil(end/2),ivar(iv));
uu = sqrt(uu) ./ utau .*sqrt(dataF(imF,1:ceil(end/2),RHO)/rho_w);
plot( yy / delBL, uu, 'k--','LineWidth',LW);
end

xlim([0 1.5])

h = ylabel('$\frac{u^\prime_{i_{rms}}}{u_{\tau}} \sqrt{\frac{\overline{\rho}}{\rho_w}}$'); set(h,'Interpreter','latex','FontSize',FSn);
h = xlabel('$y/\delta$'); set(h,'Interpreter','latex','FontSize',FSn);
%h = legend('Mesh A','Mesh B','Mesh C','Location','NorthWest');
%set(h,'Interpreter','latex','FontSize',FSn);
legend boxoff;
box on;

if (pdfE)
fig2pdf([path,f(4).name],4);
end


end
