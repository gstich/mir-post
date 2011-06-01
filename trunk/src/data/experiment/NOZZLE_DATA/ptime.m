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
    tap(i).PSD = sqrt((tap(i).p - tap(i).pbar).^2)./tap(i).pbar ;
end

figure(1);
for i=1:5
    hold all;   
    plot(tap(i).p(1000:4000))
end

figure(2)
L=10000;
Fs=srate;
hw = hann(L,'periodic');

for i=2:5

y = tap(i).PSD(1:L);
y = y.*hw;
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
%Y = fftshift(fft(y,NFFT))/L;
Y = (fft(y,NFFT))/L;
f = Fs/2*linspace(0,1,NFFT/2+1);

c(i-1,:) = 2*abs(f'.*Y(1:NFFT/2+1));
% Plot single-sided amplitude spectrum.
figure(i)
semilogx(f,2*abs(f'.*Y(1:NFFT/2+1)));hold all; 
%loglog(f,2*abs(f'.*Y(1:NFFT/2+1)));hold all; 
%xlim([f(ceil(size(f,2)/2)) f(end) ])
end
title('Single-Sided Amplitude Spectrum of y(t)')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')



%% Reduce data to a mere N points
N = 1000;
f1 = f(2);
f2 = f(end);
FF = logspace(log10(f1),log10(f2),N);
for i=1:4
    b(i,:) = interp1(f,c(i,:),FF);
    b(i,:) = b(i,:) / sum(b(i,:));
end


% Make a 2d Contour plot of exp. data
lfreq = size(FF,2);
x = zeros(4,lfreq);
y = zeros(4,lfreq);
for i=1:lfreq
    x(:,i) = [0, 1.27, 2.54, 3.81];
    % Filter this term
    %for k=1:4
    %    for ii=2:min(size(b))-1
    %        b(ii,:)=(b(ii-1,:)+b(ii+1,:))/2;
    %    end
    %end
end

for i=1:4
    y(i,:) = log10(FF);
    % Filter this term
    for k=1:2
        for ii=2:max(size(b))-1
            b(:,ii)=(b(:,ii-1)+b(:,ii+1))/2;
        end
    end
end

figure(7);
[cp,cp] = contourf(x,y,b,12);
set(cp,'edgecolor','none');
xray = flipud(gray);
colormap(xray)




