clear all;
clc;
close all;

dir = 'TIME-RESOLVED';
srate = 200e3;
V = 1;
Pamb = 14.7;

tap(1).name = [ dir , '/SOE00316.PR1'];
tap(1).m =                                  3.030500	;
tap(1).b =                                  -0.798400	;


tap(2).name = [ dir , '/SOE00316.PR2'];
tap(2).m =                                  2.948200	;
tap(2).b =                                  -0.715200	;

tap(3).name = [ dir , '/SOE00316.PR3'];
tap(3).m =                                  2.944700	;
tap(3).b =                                  0.524300	;

tap(4).name = [ dir , '/SOE00316.PR4'];
tap(4).m =                                  2.943500	;
tap(4).b =                                  -1.257100	;

tap(5).name = [ dir , '/SOE00316.PR5'];
tap(5).m =                                  2.922600	;
tap(5).b =                                  -0.596500;


for i=1:5
    tap(i).praw = load(tap(i).name);
    tap(i).p = Pamb + tap(i).praw*tap(i).m + tap(i).b;
    tap(i).pbar = mean( tap(i).p );
    tap(i).PSD = sqrt((tap(i).p - tap(i).pbar).^2);
end

figure(1);
for i=1:5
    hold all;   
    plot(tap(i).p(1000:4000))
end

figure(2)
L=1000;
Fs=srate;
hw = hann(L,'periodic');
for i=2:5

y = tap(i).PSD(1:L);
y = y.*hw;
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
%Y = fftshift(fft(y,NFFT))/L;
Y = (fft(y,NFFT))/L;
f = Fs/2*linspace(0,1,NFFT/2+1);

% Plot single-sided amplitude spectrum.
figure(i)
semilogx(f,2*abs(f'.*Y(1:NFFT/2+1)));hold all; 
%loglog(f,2*abs(Y(1:NFFT/2+1)));hold all; 
xlim([f(ceil(size(f,2)/2)) f(end) ])
end
title('Single-Sided Amplitude Spectrum of y(t)')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')


