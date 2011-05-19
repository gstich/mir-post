clear all;
clc;
close all;


%% Plot the BL data: Mean and Reynolds stress
pdfE = true;
res(1).name = 'coarse';
res(2).name = 'medium';
ref = 2;

% Figure name
figs(1).name = '2dspancorr';
figs(2).name = '1dspancorr';

% Plot symbols
res(1).sym = '--';
res(2).sym = '-.';

% Figure option
LW = 2;         % LineWidth
FSn = 18;       % FontSize labels
FSa = 12;       % FontSize axis
cmax = .2;      % Contour Level Max
Nc = 64;        % Number of contour level
Ht = 1.78;

%% Coarse
figure(1);
load('span_corr/cor.mat');
XX = XX/Ht;
YY = YY/Ht;
contourf(XX,YY+1,tay,[0:.01:.2],'edgecolor','none')


%% Correlation plots
figure(2);hold on;
x = [1:1:size(Cor,1)];

plot(  x/x(end) *.5, Cor(:,1) ,'--','LineWidth',LW);
plot(  x/x(end) *.5, Cor(:,106) ,'--','LineWidth',LW);
plot(  x/x(end) *.5, Cor(:,399), '--','LineWidth',LW);

%% Medium
figure(1);hold on;
load('span_corr/med.mat');
XX = XX/Ht;
YY = YY/Ht;
contourf(XX,YY-1,tay,[0:cmax/Nc:cmax],'edgecolor','none')
colorbar;
axis equal
box on;
h1 = xlabel(['$x/ H_t$']);
set(h1,'Interpreter','latex','FontSize',FSn);
h2 = ylabel('$ y/H_t $');
set(h2,'Interpreter','latex','FontSize',FSn);
set(gca,'FontSize',FSa);



%% Correlation plots
figure(2);hold on;
x = [1:1:size(Cor,1)];
plot(  x/x(end) *.5, Cor(:,1) ,'LineWidth',LW);
plot(  x/x(end) *.5, Cor(:,106) ,'LineWidth',LW);
plot(  x/x(end) *.5, Cor(:,237) ,'LineWidth',LW);
box on;
h1 = xlabel(['$z/ H_t$']);
set(h1,'Interpreter','latex','FontSize',FSn);
h2 = ylabel('$ R_{11}(u) $');
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



