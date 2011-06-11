function  spectral


%res = 'medium';
res = 'coarse';


switch(res)
    case('coarse');
        Bfile = 'DATA/ptot_plume.cor';
        npts = 19;
        ofile = '../data/spectra/separation/spectra_coarse.mat';
    case('medium')    
        Bfile = 'DATA/ptot_plume.med';
        npts = 19;
        ofile = '../data/spectra/separation/spectra_medium.mat'   ;     
end


b = 0;
for i=1:npts
    file = [Bfile, int2str(i) ];
    data = load(file);
    time = data(:,1);
    p = data(:,5);
    u = data(:,2);
    v = data(:,3);
    w = data(:,4);
    rho = data(:,7);
    ptot = p + rho/2.*(u.*u + v.*v + w.*w);
    
    %ptot = 2.2*ptot;
    a = get_spectra(ptot);
    b = b + a;
end

b = b/npts;


tt = max(size(p))*5e-6 * 10;
c = b(2,:).^2*tt;
b(2,:) = c ;%/ max(c);
loglog(b(1,:),b(2,:),'k');hold all; 

psd_e = load('../data/experiment/psd_exp.dat');
%psd_e(:,2) = psd_e(:,2)/psd_e(1,2);
loglog(psd_e(:,1),psd_e(:,2),'ro');

xlabel('f H_t/U_p')
ylabel('Pa^2/Hz')

a = 1*10^3.5;
m = -25/9;
x = logspace(-1,.5,50);
y = a*x.^m;

%hold on;
%plot(x,y,'k--');


save(ofile, 'b', 'psd_e');

%% Compare to exp. data
%plot_exp;
%compare_exp(x,y,'k--');
%compare_exp(b(1,:),2*b(2,:).^2*tt,'b--');



end


function compare_exp(x,y,sym)
a = 1;
% Read in the exp. data file

% Data extents on plot
xx = [-2.5,0];
yy = [2,6];
% Data extents on picture-in pixels
LL = [146,493];
UR = [716, 57];


xdata = log10(x);
ydata = log10(y);

xdata = (xdata - xx(1) ) / (xx(2)-xx(1));
xdata = xdata * (UR(1)-LL(1)) + LL(1);

ydata = -(ydata - yy(2)) / (yy(2)-yy(1));
ydata = ydata * (LL(2)-UR(2)) + UR(2);

hold on;
plot(xdata,ydata,sym,'LineWidth',2);
end

function plot_exp
a = 1;
figure(3);
exp = imread('ptot_spectra.png');

image(exp);
end


function a = get_spectra(ptot)

Ht = 1.78;
Up = 32940*1.603;
L=size(ptot,1);
Fs = 200e3;
hw = hanning(L,'periodic');

y = ptot - mean(ptot);
y = y/10;   % Convert to Pascals
%y = y.*y;
y = y.*hw;


NFFT = 2^nextpow2(L); % Next power of 2 from length of y
figure(2)
Y = (fft(y,NFFT))/L;
f = Fs/2*linspace(0,1,NFFT/2+1);

%% Divide by frequency

a(1,:) = f*Ht/Up;
a(2,:) = 2*abs(Y(1:NFFT/2+1));  %./f';




end