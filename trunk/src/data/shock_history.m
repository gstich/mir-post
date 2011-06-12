clear all;
clc;
%close all;

% Figure name
figs(1).name = 'shock_history';
figs(2).name = 'shock_spectra';
figs(3).name = 'shock_compare';

pdfE = false;

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
%plot(timeC(1:f-s+1),Xs_cor(s:f,1),timeC(1:f-s+1),Xs_cor(s:f,2),timeC(1:f-s+1),Xs_cor(s:f,3));
%plot(timeM,Xs_med(:,1),timeM,Xs_med(:,2),timeM,Xs_med(:,3));

%plot(timeC,Xs_cor(:,1),'LineWidth',LW);
%plot(timeM,Xs_med(:,1),'k','LineWidth',LW);
plot(time*Up/Ht/1e3,XSS(:,1),'k','LineWidth',LW);hold on;
plot(time*Up/Ht/1e3,XSS(:,2),'b','LineWidth',LW);hold on;
plot(time*Up/Ht/1e3,XSS(:,3),'r','LineWidth',LW);

%figure(3);
%plot(time,abs(XSS(:,3)-XSS(:,1)),'k');
%xlim([0 .6]);
ylim([-.6 .6]);
box on;
h1 = xlabel(['$tU_p/H_t$']);
set(h1,'Interpreter','latex','FontSize',FSn);
h2 = ylabel('$X_s / H_t$');
set(h2,'Interpreter','latex','FontSize',FSn);
set(gca,'FontSize',FSa);
h3 = legend('Bottom','Top','Centerline');
set(h3,'Interpreter','latex','FontSize',FSn);
legend boxoff;


%shock = Xs_cor(:,1);
%shock = Xs_med(:,1);
shock = XSS(:,1);


L=size(shock,1);
Fs = 200e3;
hw = hanning(L,'periodic');





NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = zeros(NFFT,1);
figure(2);
for i=1:3
    %shock = Xs_med(:,i);
    %shock = Xs_cor(:,i);
    shock = XSS(:,i);
    y = shock - mean(shock);
    y = y/100; % Convert to meters
    y = y.*hw;
    Y = Y + (fft(y,NFFT))/L;
    
end
f = Fs/2*linspace(0,1,NFFT/2+1);

%% Divide by frequency
a(1,:) = f*Ht/Up;
a(2,:) = 2*abs(Y(1:NFFT/2+1)) / 3;

figure(2);
loglog( a(1,:), a(2,:).*a(2,:) ,'r-.','LineWidth',LW)

hold on;
eoff = 1;
psd_e = load('experiment/shock_psd.dat');
%psd_e(:,2) = psd_e(:,2)/psd_e(1,2);
loglog(psd_e(eoff:end,1),psd_e(eoff:end,2),'k--','LineWidth',LW);


% Read in the exp. data file		
			
% Data extents on plot		
xx = [0,20];		
yy = [-0.5,0.7];		
		
xoff = .55;

% Data extents on picture-in pixels		
LL = [126,464];%[509,348];		
UR = [1479,23];%[712, 14];		
		
figure(3);		
exp = imread('shock_history/shock_history.png');		
		
image(exp);		
xpix = size(exp,2);		
ypix = size(exp,1);		
%		
%		
xdata = time;		
ydata = XSS(:,1);		
		
xdata = (xdata - xx(1) ) / (xx(2)-xx(1));		
xdata = xdata * (UR(1)-LL(1)) + LL(1);		
%		
ydata = -(ydata - yy(2)) / (yy(2)-yy(1));		
ydata = ydata * (LL(2)-UR(2)) + UR(2);		
%		
pxoff = xoff*xpix;
hold on;		
plot(xdata+pxoff,ydata,'b-','LineWidth',LW);

axis off;


figure(2);
xlim([.005 2]);
ylim([10^-12 10^-4.5]);
box on;
h1 = xlabel(['$tU_p/H_t$']);
set(h1,'Interpreter','latex','FontSize',FSn);
h2 = ylabel('$S_{xx}(m^2/Hz)$');
set(h2,'Interpreter','latex','FontSize',FSn);
set(gca,'FontSize',FSa);


%a = 1;
%b = 1;
%c = xpix;
%d = ypix*.8;
%set(gcf,'Position',[a b c d])
%set(gcf, 'PaperSize', [10.5 3]);


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
