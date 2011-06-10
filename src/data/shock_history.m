clear all;
clc;
%close all;

% Figure name
figs(1).name = 'shock_history';

pdfE = false;

res = 'coarse';


% Figure option
LW = 2;         % LineWidth
FSn = 18;       % FontSize labels
FSa = 12;       % FontSize axis


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

XSS = (-11.7 + XSS) / 1.78 + 3.5;
%Xs_cor = (-11.7 + Xs_cor) / 1.78 + 3.5;
%Xs_med = (-11.7 + Xs_med) / 1.78 + 3.5;

%steps = size(Xs_cor,1);
%timeC = linspace(0,(steps-1)*dt,steps);
%timeC = timeC * 1000;

%steps = size(Xs_med,1);
%timeM = linspace(0,(steps-1)*dt,steps);
%timeM = timeM * 1000;

steps = size(XSS,1);
time = linspace(0,(steps-1)*dt,steps);
time = time * 1000;


s = 200;
f = size(Xs_cor,1);

figure(1);hold on;
%plot(timeC(1:f-s+1),Xs_cor(s:f,1),timeC(1:f-s+1),Xs_cor(s:f,2),timeC(1:f-s+1),Xs_cor(s:f,3));
%plot(timeM,Xs_med(:,1),timeM,Xs_med(:,2),timeM,Xs_med(:,3));

%plot(timeC,Xs_cor(:,1),'LineWidth',LW);
%plot(timeM,Xs_med(:,1),'k','LineWidth',LW);
plot(time,XSS(:,1),'k','LineWidth',LW);hold on;
plot(time,XSS(:,2),'b','LineWidth',LW);hold on;
plot(time,XSS(:,3),'r','LineWidth',LW);

%figure(3);
%plot(time,abs(XSS(:,3)-XSS(:,1)),'k');

box on;
h1 = xlabel(['Time (msec)']);
set(h1,'Interpreter','latex','FontSize',FSn);
h2 = ylabel('$X_s / H_t$');
set(h2,'Interpreter','latex','FontSize',FSn);
set(gca,'FontSize',FSa);


%shock = Xs_cor(:,1);
%shock = Xs_med(:,1);
shock = XSS(:,1);

Ht = 1.78;
Up = 32940*1.93;
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
loglog( a(1,:), a(2,:).*a(2,:) ,'k')

hold on;
psd_e = load('experiment/shock_psd.dat');
%psd_e(:,2) = psd_e(:,2)/psd_e(1,2);
loglog(psd_e(:,1),psd_e(:,2),'ro');




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
