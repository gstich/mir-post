function  spectral


npts = 5;
Bfile = '../ptot_time.med';

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

    a = get_spectra(ptot);
    b = b + a;
end

b = b/npts;



loglog(b(1,:),b(2,:));hold all; 



xlabel('f H_t/U_p')
ylabel('Pa^2/Hz')

a = 3*10^2;
m = -8/3;
x = logspace(-1,.5,50);
y = a*x.^m;

hold on;
plot(x,y,'k--');

end





function a = get_spectra(ptot)

Ht = 1.78;
Up = 32940*1.93;
L=size(ptot,1);
Fs = 200e3;
hw = hanning(L,'periodic');

y = ptot - mean(ptot);
y = y/10;   % Convert to Pascals
y = y.*y;
y = y.*hw;


NFFT = 2^nextpow2(L); % Next power of 2 from length of y
figure(2)
Y = (fft(y,NFFT))/L;
f = Fs/2*linspace(0,1,NFFT/2+1);

%% Divide by frequency

a(1,:) = f*Ht/Up;
a(2,:) = 2*abs(Y(1:NFFT/2+1))./f';




end