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

Xsamp = 2.0;

pdfE = false;
pdfE = true;

plt.FSn = FSn;
plt.FSa = FSa;

plt.pdf = pdfE;


%path = '~/Dropbox/Britton/THESIS/Figures/nozzle/convergence/';
path = '/Users/bolson/Dropbox/Britton/THESIS/Figures/nozzle/convergence/';

% Figure 1- Coeff. of friction
pfig(1) = 1;
f(1).name = 'Cf_muComp';

% Figure 4- Boundary layer profile, raw
pfig(4) = 0;
f(4).name = 'UBLmean';

% Figure 5- Boundary layer profile, wall units
pfig(5) = 1;
f(5).name = 'BLlog_muComp';

% Figure 6- Artificial terms
pfig(6) = 1;
f(6).name = 'artMU_muComp';

% Figure 11- BL RMS profiles
pfig(11) = 1;
f(11).name = 'BLrms_muComp';



var_map

% Load the 2 mesh sizes
%load('../../data/jcp/fin_muScalar.mat');
load('../../data/jcp/ufinScalar.mat');
dataM = data;

%load('../../data/jcp/fin_muNew.mat');
load('../../data/jcp/ufinVector.mat');
dataF = data;


%load('../../data/jcp/fin_nomu.mat');
load('../../data/jcp/ufinNomu.mat');
dataC = data;

Ht = 0.1;
Up = data(1,end,U);
Pamb = 1e6;
delBL = 0.1;

% Get some variables

um = sqrt( dataM(:,:,U).^2 + dataM(:,:,V).^2 );
uf = sqrt( dataF(:,:,U).^2 + dataF(:,:,V).^2 );
uc = sqrt( dataC(:,:,U).^2 + dataC(:,:,V).^2 );

xm = dataM(:,:,X);ym = dataM(:,:,Y);
xf = dataF(:,:,X);yf = dataF(:,:,Y);
xc = dataC(:,:,X);yc = dataC(:,:,Y);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 1- CF convergence plot
if (pfig(1) == 1);


    % Scaling
    mu_w = dataM(:,1,MU);
    rho_w = dataM(:,1,RHO);
    rho_in = dataM(1,end,RHO);
    dudy = ( um(:,2) - um(:,1) )./ ( ym(:,2) - ym(:,1) );
    tauw = mu_w .* dudy;
    Cf = 2*tauw/(Up^2*rho_in);
    
    figure(1);
    plot(xm(:,1)/delBL,Cf,'k-.','LineWidth',LW);
    
    % Scaling
    mu_w = dataF(:,1,MU);
    rho_w = dataF(:,1,RHO);
    rho_in = dataF(1,end,RHO);
    dudy = ( uf(:,2) - uf(:,1) )./ ( yf(:,2) - yf(:,1) );
    tauw = mu_w .* dudy;
    Cf = 2*tauw/(Up^2*rho_in);
    
    hold on;
    plot(xf(:,1)/delBL,Cf,'k--','LineWidth',LW);
    
    % Scaling
    mu_w = dataC(:,1,MU);
    rho_w = dataC(:,1,RHO);
    rho_in = dataC(1,end,RHO);
    dudy = ( uc(:,2) - uc(:,1) )./ ( yc(:,2) - yc(:,1) );
    tauw = mu_w .* dudy;
    Cf = 2*tauw/(Up^2*rho_in);
    
    hold on;
    plot(xf(:,1)/delBL,Cf,'k-','LineWidth',LW);
    
    
    ylim([0 5e-3])
    
    
    box on;
    
    plt.fig = 1;
    plt.file = f(1).name;
    plt.xlabel = '$x/\delta$';
    plt.ylabel = '$C_f$';
    plt.legend = {'Medium','Fine'};
    plt.legend = {'Scalar','Directional','$C_\mu = 0$'};
    plt.misc = 'legend boxoff;';
    pretty_plot(plt)
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 4- BL profile at the in nozzle inlet
if (pfig(4) == 1);

[y,imM] = min( abs( xm(:,1)-Xsamp ));
[y,imF] = min( abs( xf(:,1)-Xsamp ));
    

figure(4);
yy = ym(imM,1:ceil(end/2));
yy = yy - yy(1);
uu = um(imM,1:ceil(end/2));hold all
%semilogx( yy / Ht, uu / Up, 'b');
plot( yy / Ht, uu / Up, 'k-.','LineWidth',LW);

yy = yf(imF,1:ceil(end/2));
yy = yy - yy(1);
uu = uf(imF,1:ceil(end/2));
%semilogx( yy / Ht, uu / Up, 'b');
plot( yy / Ht, uu / Up, 'k-','LineWidth',LW);

xlim([0 1.15])

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

[y,imM] = min( abs( xm(:,1)-Xsamp ));
[y,imF] = min( abs( xf(:,1)-Xsamp ));
[y,imC] = min( abs( xc(:,1)-Xsamp )); 

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
[a,b] = min(abs(uu-uu(end)*.99));
delBLm = yy(b)

ut1 = utau;
tw1 = tauw;

% Van Driest Transform
uvd = zeros(1,max(size(uu)));
for j=1:max(size(uu))
  uvd(j) = trapz( uu(1:j),sqrt( dataM(imM,1:j,RHO)/rho_w ) );
end 



semilogx( yy / del, uvd / utau, 'k-.','LineWidth',LW);hold all;

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
[a,b] = min(abs(uu-uu(end)*.99));
delBLf = yy(b)

abs(tw1-tauw)/tauw
abs(ut1-utau)/utau

% Van Driest Transform
uvd = zeros(1,max(size(uu)));
for j=1:max(size(uu))
  uvd(j) = trapz( uu(1:j),sqrt( dataF(imF,1:j,RHO)/rho_w ) );
end 

semilogx( yy / del, uvd / utau, 'k--','LineWidth',LW);



yy = yc(imC,1:ceil(end/2));
yy = yy - yy(1);
% Scaling
mu_w = dataC(imC,1,MU);
rho_w = dataC(imC,1,RHO);
dudy = ( uc(imC,2) - uc(imC,1) )/ ( yc(imC,2) - yc(imC,1) );
tauw = mu_w * dudy
utau = sqrt( tauw / rho_w );
del = mu_w / ( rho_w * utau );
uu = uc(imC,1:ceil(end/2));
[a,b] = min(abs(uu-uu(end)*.99));
delBLc = yy(b)

abs(tw1-tauw)/tauw
abs(ut1-utau)/utau

% Van Driest Transform
uvd = zeros(1,max(size(uu)));
for j=1:max(size(uu))
  uvd(j) = trapz( uu(1:j),sqrt( dataC(imC,1:j,RHO)/rho_w ) );
end 

semilogx( yy / del, uvd / utau, 'k-','LineWidth',LW);



x1 = logspace(-.5,1.2,20);
x2 = logspace(.8,3,20);
k = .41;
C = 5.1;
hold on;
semilogx(x1,x1,'k--');hold on;
semilogx(x2,1/k*log(x2)+C,'k--');

xlim([.5 700])
plt.fig = 5;
plt.ylabel = '$u_{vd}^+$';
plt.xlabel = '$y^+$';
%plt.legend = {'$C_\mu = .002$','$C_\mu = 0$','Location','NorthWest'};
plt.legend = {'Scalar','Directional','$C_\mu = 0$','Location','NorthWest'};
%plt.legend = {'Medium','Fine','Location','NorthWest'};
plt.misc = 'legend boxoff;';
plt.file = f(5).name;

pretty_plot(plt);


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 6- 
if (pfig(6) == 1);
figure(6);
    
[y,imM] = min( abs( xm(:,1)-Xsamp ));
[y,imF] = min( abs( xf(:,1)-Xsamp ));
%[y,imC] = min( abs( xc(:,1)-Xsamp ));

vf2 = dataF(imF,:,MU);
yif = yf(imF,:);
vm2 = dataM(imM,:,MU);
yim = ym(imM,:);
%vc2 = dataC(imC,:,MU);
%yic = yc(imC,:);


vm = dataM(imM,:,MUc);

%Scalar
loglog(yim/Ht,vm./vm2-1,'k-','Linewidth',LW);hold on;


% Directional
vf = dataF(imF,:,MUa);
loglog(yif/Ht,vf./vf2-1,'k--','Linewidth',LW);

vf = dataF(imF,:,MUb);
loglog(yif/Ht,vf./vf2-1,'k-.','Linewidth',LW);

vf = dataF(imF,:,MUc);
loglog(yif/Ht,vf./vf2-1,'k:','Linewidth',LW);


vf = dataF(imF,:,MUa);
vm = dataM(imM,:,MUa);

%plot(yim/Ht,vm./vm2,'b--','Linewidth',LW/2);
%plot(yif/Ht,vf./vf2,'r--','Linewidth',LW);
xlim([0 2.5])
ylim([1e-5 100])

%h = legend('$\mu^*_{\eta}$ Mesh A','$\mu^*_{\eta}$ Mesh B');
%set(h,'Interpreter','latex','FontSize',FSn);

plt.fig = 6;
plt.xlabel = '$y/\delta$';
plt.ylabel = '$\mu*/\mu$';
plt.legend = {};
plt.legend = {'$\mu^*$','$\mu^*_{x}$','$\mu^*_{y}$','$\mu^*_{z}$'};
plt.misc = 'legend boxoff; box on;';
plt.file = f(6).name;

pretty_plot(plt);


end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 11- BL profile at the in nozzle inlet
if (pfig(11) == 1);

[y,imM] = min( abs( xm(:,1)-Xsamp ));
[y,imF] = min( abs( xf(:,1)-Xsamp ));
    
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
yy = ym(imM,:);
yy = yy - yy(1);
uu = dataM(imM,:,ivar(iv));
uu = sqrt(uu) ./ utau .*sqrt(dataM(imM,:,RHO)/rho_w);
hold all;
plot( yy / delBLm, uu, 'k-.','LineWidth',LW);
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
yy = yf(imF,:);
yy = yy - yy(1);
uu = dataF(imF,:,ivar(iv));
uu = sqrt(uu) ./ utau .*sqrt(dataF(imF,:,RHO)/rho_w);
plot( yy / delBLf, uu, 'k--','LineWidth',LW);
end

% Scaling
mu_w = dataC(imC,1,MU);
rho_w = dataC(imC,1,RHO);
dudy = ( uc(imC,2) - uc(imC,1) )/ ( yc(imC,2) - yc(imC,1) );
tauw = mu_w * dudy
utau = sqrt( tauw / rho_w );
del = mu_w / ( rho_w * utau );
cf = 2*tauw / dataC(imC,end,RHO)/uc(imC,end)^2

for iv = 1:max(size(ivar));
yy = yc(imC,:);
yy = yy - yy(1);
uu = dataC(imC,:,ivar(iv));
uu = sqrt(uu) ./ utau .*sqrt(dataC(imC,:,RHO)/rho_w);
plot( yy / delBLc, uu, 'k-','LineWidth',LW);
end


xlim([0 1.5])

h = ylabel('$\frac{u^\prime_{i_{rms}}}{u_{\tau}} \sqrt{\frac{\overline{\rho}}{\rho_w}}$'); set(h,'Interpreter','latex','FontSize',FSn);
h = xlabel('$y/\delta$'); set(h,'Interpreter','latex','FontSize',FSn);
%h = legend('Mesh A','Mesh B','Mesh C','Location','NorthWest');
%set(h,'Interpreter','latex','FontSize',FSn);
legend boxoff;
box on;

plt.fig = 11;
%plt.pdf = true;
plt.ylabel = '$\frac{u^\prime_{i_{rms}}}{u_{\tau}} \sqrt{\frac{\overline{\rho}}{\rho_w}}$';
plt.xlabel = '$y/\delta$';
plt.legend = {};
%plt.legend = {'Scalar','Directional','No SGS','Location','NorthEast'};
plt.misc = 'legend boxoff; box on;';
plt.file = f(11).name;

pretty_plot(plt);



end





