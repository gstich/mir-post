clear all;
clc;
close all;

% Figure name
figs(1).name = 'shock_history';

pdfE = false;
% Figure option
LW = 2;         % LineWidth
FSn = 18;       % FontSize labels
FSa = 12;       % FontSize axis


pref = 7.5e5;
dt = 5.0e-6;
Ht = 1.78;

%% Load the data
load shock_history/Xs_cor;
Xs_cor = XSS;
%Xs_cor = Xs_cor / Ht + 3;

load shock_history/Xs_med;
Xs_med = XSS;
%Xs_med = Xs_med / Ht + 3;

Xs_cor = (-11.7 + Xs_cor) / 1.78 + 3.5;
Xs_med = (-11.7 + Xs_med) / 1.78 + 3.5;

steps = size(Xs_cor,1);
timeC = linspace(0,(steps-1)*dt,steps);
timeC = timeC * 1000;

steps = size(Xs_med,1);
timeM = linspace(0,(steps-1)*dt,steps);
timeM = timeM * 1000;


s = 200;
f = size(Xs_cor,1);

figure(1);hold on;
%plot(timeC(1:f-s+1),Xs_cor(s:f,1),timeC(1:f-s+1),Xs_cor(s:f,2),timeC(1:f-s+1),Xs_cor(s:f,3));
%plot(timeM,Xs_med(:,1),timeM,Xs_med(:,2),timeM,Xs_med(:,3));

plot(timeC(1:f-s+1),Xs_cor(s:f,1),'LineWidth',LW);
plot(timeM,Xs_med(:,1),'k','LineWidth',LW);

box on;
h1 = xlabel(['Time (msec)']);
set(h1,'Interpreter','latex','FontSize',FSn);
h2 = ylabel('$X_s / H_t$');
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
