function ff = compute_spectra(F)

ptot = F;


L=size(ptot,1);
Fs = 200e3;
hw = hanning(L,'periodic');

y = ptot - mean(ptot);
y = y.*y;
y = y.*hw;


NFFT = 2^nextpow2(L) % Next power of 2 from length of y
Y = (fft(y,NFFT))/L;
f = 1/2*linspace(0,1,NFFT/2+1);

%% Divide by frequency

a(1,:) = f;%*Ht/Up;
a(2,:) = 2*abs(Y(1:NFFT/2+1));  %./f';

ff = a;



end