clear all;
clc;
close all;

%% Plot the BL data: Mean and Reynolds stress
pdfE = true;
res(1).name = 'coarse';
res(1).xoff = -.24;
res(2).name = 'medium';
res(2).xoff = -.24;
ref = 2;

% Figure name
figs(1).name = 'Pcenter';
figs(2).name = 'Pwall';
figs(3).name = 'dPdx';

% Plot symbols
res(1).sym = '--';
res(2).sym = '-.';

% Figure option
LW = 2;         % LineWidth
FSn = 18;       % FontSize labels
FSa = 12;       % FontSize axis

% Data format key
key.xcen = 1;       key.pcen = 2;        key.xwall = 3;
key.pwalld = 4;     key.pwallu = 5;      
key.pwall(1) = 5;
key.pwall(2) = 4;


%% Experimental
efile = '../data/experiment/NOZZLE_DATA/MEAN/CENTERLINE/MeanCentPress.txt';
PexpC = load(efile);
figure(1);hold on;
plot(PexpC(:,1),PexpC(:,2),'r-','LineWidth',LW)

efile = '../data/experiment/NOZZLE_DATA/MEAN/WALL/MeanWallPress.txt';
PexpW = load(efile);
figure(2);hold on;
plot(PexpW(:,1),PexpW(:,2),'ro','LineWidth',LW)



% Load the data
for i=1:size(res,2)
    file = ['pressure/',res(i).name, '.dat'];
    res(i).data = load(file);
end


%% Mean Profiles
for i=1:size(res,2)
    x = res(i).data(15:end,key.xcen) + res(i).xoff;
    p = res(i).data(15:end,key.pcen);
    figure(1);hold on;
    plot(x,p,res(i).sym,'LineWidth',LW);hold on;
    
    figure(2);hold on;
    x = res(i).data(15:end,key.xwall)+ res(i).xoff;
    
    %if(i==2)
    p = res(i).data(15:end,key.pwall(i));
    %else
    %p = res(i).data(15:end,key.pwallu);
    %end
    plot(x,p,res(i).sym,'LineWidth',LW);hold on;

end

figure(1);
xlim([0 8]);
box on;
h1 = xlabel(['$x/H_t$']);
set(h1,'Interpreter','latex','FontSize',FSn);
h2 = ylabel('$P/P_0$');
set(h2,'Interpreter','latex','FontSize',FSn);
set(gca,'FontSize',FSa);


figure(2);
xlim([0 6]);
box on;
h1 = xlabel(['$x/H_t$']);
set(h1,'Interpreter','latex','FontSize',FSn);
h2 = ylabel('$P/P_0$');
set(h2,'Interpreter','latex','FontSize',FSn);
set(gca,'FontSize',FSa);


%% Plot dP/dx from x/H = [-1,2]
P0 = 1e6;
tauw = 2388.4;
%del = .0358/2;
del = .0213;
H = 1.78;
beta = P0*del/(tauw*H);
figure(3);
for i=1:size(res,2)
    %x = res(i).data(15:end,key.xcen) + res(i).xoff;
    %p = res(i).data(15:end,key.pcen);
    x = res(i).data(15:end,key.xwall) + res(i).xoff;
    p = res(i).data(15:end,key.pwall(i));
    
    dPdx = ddx(p)./ddx(x);
    plot(x,dPdx*beta,res(i).sym,'LineWidth',LW);hold on;
end
    
%p = PexpC(:,2);
%x = PexpC(:,1);
p = PexpW(:,2);
x = PexpW(:,1);

dPdx = ddx(p)./ddx(x);
%for i=1:100;dPdx = gfilter(dPdx);end;
plot(x,dPdx*beta,'ro','LineWidth',LW);hold on;

% Spalart DNS
x = linspace(-1,3,100);
p = x*0 + -.3;
plot(x,p,'k--');

figure(3);
xlim([-1,2.5])
box on;
h1 = xlabel(['$x/H_t$']);
set(h1,'Interpreter','latex','FontSize',FSn);
h2 = ylabel('$\frac{dP}{dx} \frac{\delta^*}{\tau_w}$');
set(h2,'Interpreter','latex','FontSize',FSn);
set(gca,'FontSize',FSa);



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

