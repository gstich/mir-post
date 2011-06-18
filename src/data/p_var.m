clear all;
clc;
close all;

% Figure name
figs(1).name = 'p_var';

pdfE = false;
pdfE = true;

res = 'coarse';
%res = 'medium';


% Figure option
LW = 2;         % LineWidth
FSn = 25;       % FontSize labels
FSa = 18;       % FontSize axis


pref = 7.5e5;
dt = 5.0e-6;
Ht = 1.78;

%% Load the data
load('experiment/NOZZLE_DATA/var_comp.mat');

figure(1);hold on;
plot(sX,var,'k','LineWidth',LW);
plot(eX,evar(2:5),'ko','LineWidth',LW);


xlim([2 6]);
ylim([0 .075]);
box on;
h1 = xlabel(['$X/H_t$']);
set(h1,'Interpreter','latex','FontSize',FSn);
h2 = ylabel('$\langle p^\prime p^\prime\rangle/\langle p\rangle^2$');
set(h2,'Interpreter','latex','FontSize',FSn);
set(gca,'FontSize',FSa);
legend boxoff;

h = figure(1);
%set(h,'Position',[1000,395,560*1.5,420*1.5])
%set(h,'PaperSize',[8,7]);
%set(h,'PaperPositionMode','auto');



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
