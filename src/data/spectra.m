clear all;
clc;
close all;

%% Plot the BL data: Mean and Reynolds stress
pdfE = true;
%pdfE = false;
res(1).name = 'coarse';
res(2).name = 'medium';
res(3).name = 'fine';
ref(1) = 1;
ref(2) = 2;
ref(3) = 2;

% Figure name
figs(1).name = 'PSD_compare';

% Plot symbols
res(1).sym = 'k-';
res(2).sym = 'k-.';
res(3).sym = 'k-';

% Figure option
LW = 2;         % LineWidth
FSn = 25;       % FontSize labels
FSa = 18;       % FontSize axis

% Load the data
%for i=1:size(res,2)
for i=1:1
    file = ['spectra/separation/spectra_',res(i).name, '.mat'];
    load(file);
    %res(i).data = b;
    %exp = psd_e;
    
    loglog(b(1,:),b(2,:),res(i).sym,'LineWidth',LW);hold on;
end

loglog(psd_e(:,1),psd_e(:,2),'ko','LineWidth',LW);hold on;



figure(1);

h1 = xlabel(['$fH_t/U_p$']);
set(h1,'Interpreter','latex','FontSize',FSn);
h2 = ylabel('$PSD^*$');
set(h2,'Interpreter','latex','FontSize',FSn);
box on;
set(gca,'FontSize',FSa);
xlim([10^-3 3]);
ylim([10^-6 2]);
set(gca,'Position',[.13,.13,.775,.8107])


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

