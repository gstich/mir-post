clear all;
clc;
close all;


%% Plot the BL data: Mean and Reynolds stress
pdfE = false; %
pdfE = true;
res(1).name = 'coarse';
res(2).name = 'medium';
res(3).name = 'fine';

% Plot symbols
res(1).Lz = 128*.021;
res(2).Lz = 256*.014;
res(3).Lz = 384*.0105;

% Figure name
figs(1).name = 'spec_span_u';
figs(2).name = 'corr_span_p';

% Plot symbols
res(1).sym = '--';
res(2).sym = '-.';
res(3).sym = '-';

% Figure option
LW = 2;         % LineWidth
FSn = 25;       % FontSize labels
FSa = 18;       % FontSize axis
cmax = .2;      % Contour Level Max
Nc = 64;        % Number of contour level
Ht = 1.78;


%% Plot 1- u' spectra in LB
var = 'u';
for i=1:3
    fdat = ['span_spec/',res(i).name,'-',var,'.mat'];
    load(fdat);
    res(i).spec = spec;
    res(i).corr = corr;
end

%% spec ( k , npts, 2)

figure(1);
for i=1:3
    nn = size(res(i).spec,1)-1;
    knorm = linspace(1,nn,nn)/res(i).Lz;
    loglog(knorm,res(i).spec(2:end,4,2),['k',res(i).sym],'LineWidth',LW);hold on;
end

figure(1);
for i=1:3
    nn = size(res(i).spec,1)-1;
    knorm = linspace(1,nn,nn)/res(i).Lz;
    loglog(knorm,res(i).spec(2:end,3,2),['b',res(i).sym],'LineWidth',LW);hold on;
end



xlim([.25 65])
ylim([10^5 2*10^8])
box on;
h1 = xlabel(['$k_z$']);
set(h1,'Interpreter','latex','FontSize',FSn);
h2 = ylabel('$ E_u(k) $');
set(h2,'Interpreter','latex','FontSize',FSn);
set(gca,'FontSize',FSa);

h = figure(1);
set(h,'Position',[300,395,560*1.5,420*1.5])
set(h,'PaperSize',[8,7]);
set(h,'PaperPositionMode','auto');


%% Figure 2: Corr of pressure
%% Plot 1- u' spectra in LB
var = 'p';
for i=1:3
    fdat = ['span_spec/',res(i).name,'-',var,'.mat'];
    load(fdat);
    res(i).spec = spec;
    res(i).corr = corr;
end

% Corr( Z, npts)
figure(2);hold on;
for i=1:3
    nn = size(res(i).corr,1) / 2;
    knorm = linspace(0,.5,nn);
    plot(knorm,res(i).corr(1:end/2,4),['k',res(i).sym],'LineWidth',LW);
end

for i=1:3
    nn = size(res(i).corr,1) / 2;
    knorm = linspace(0,.5,nn);
    plot(knorm,res(i).corr(1:end/2,3),['b',res(i).sym],'LineWidth',LW);
end

plot(knorm,knorm*0,'k--','LineWidth',1);

xlim([0 .5])
ylim([-.35 1])
box on;
h1 = xlabel(['$z/L_z$']);
set(h1,'Interpreter','latex','FontSize',FSn);
h2 = ylabel('$ R_{11}(p) $');
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



