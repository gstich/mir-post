clear all;
clc;
%close all;

% Figure name
figs(1).name = 'mech_phase';

pdfE = false;

%res = 'coarse';
res = 'medium';


% Figure option
LW = 2;         % LineWidth
FSn = 25;       % FontSize labels
FSa = 18;       % FontSize axis


pref = 7.5e5;
dt = 5.0e-6;
Ht = 1.78;

%% Load the data
switch(res)
    case('coarse')
        load shock_history/coarse.mat;
        Xs_cor = XSS;
        %Xs_cor = Xs_cor / Ht + 3;
    case('medium')
        load shock_history/medium.mat;
        Xs_med = XSS;
        %Xs_med = Xs_med / Ht + 3;
end

XSS = (-11.7 + XSS) / 1.78 + 3.55;


steps = size(XSS,1);
time = linspace(0,(steps-1)*dt,steps);
time = time * 1000;

Ht = 1.78;
Up = 32940*1.603;

s = 200;
f = size(XSS,1);

figure(1);hold on;
toff = 550;

t1 = toff;
tf = 1213;
time = time*Up/Ht/1e3;
time2 = time(t1:tf)-time(t1);
shk = XSS(t1:tf,3);

for i=1:50
shk = gfilter(shk);
end

plot(time2,shk,'k','LineWidth',LW);hold on;


% Load the exit plane data

dt = time2(2)-time2(1);
exit_data = load('mechanism/exit_plane.med.full');
exit_data = load('mechanism/exit_plane.med.full.wide');
time3 = exit_data(:,1)-exit_data(1,1) - toff;
time3 = time3*dt;

pexit = exit_data(:,5);
pave = mean(pexit);
pexit = (pexit-pave)/pave;
for i=1:50
pexit = gfilter(pexit);
end

uexit = exit_data(:,2);
uave = mean(uexit);
uexit = (uexit-uave)/uave;
for i=1:50
uexit = gfilter(uexit);
end

figure(1);
plot(time3,pexit*10,'b','LineWidth',LW);hold on;
plot(time3,uexit*5,'r','LineWidth',LW);hold on;


%% Load the separation region 
toff = 0;
sep_data = load('mechanism/sep_length.dat');
time4 = (sep_data(:,1)-sep_data(1,1) + toff)*dt;
delU = sep_data(:,2);
delL = sep_data(:,3);
core = delU-delL;
core = (core-mean(core)) * 2 ;%/(1.78*1.6); 
core = gfilter(core);

figure(1);hold on;
plot(time4,core,'g','LineWidth',LW);

xline = linspace(0,100,10);
figure(1);hold on;
plot(xline,0*xline,'k-.','LineWidth',1)

xlim([0 66]);
ylim([-1. 1.]);
box on;
h1 = xlabel(['$tU_p/H_t$']);
set(h1,'Interpreter','latex','FontSize',FSn);
h2 = ylabel('$X_s / H_t$');
set(h2,'Interpreter','latex','FontSize',FSn);
set(gca,'FontSize',FSa);
h3 = legend('Shock','$P_{exit}$','$U_{exit}$','$A_{eff}$');
set(h3,'Interpreter','latex','FontSize',FSn);
legend boxoff;



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
