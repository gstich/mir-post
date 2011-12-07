clear all;
clc;
close all;

path(path,'../../../matlab_scripts');
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
    plot(tap(i).p(1000:20000))
end


figure(9);
t1=1;
tf=20000;
xp = [0, 1.27, 2.54, 3.81] + 2.42;  % 1st and last probes are placed .5in apart (1.27cm)
xp = [0, 0.4233, .8467, 1.27] + 2.42;

ps1 = mean(tap(2).p);  % Mean left bound pressure
ps2 = mean(tap(5).p);  % Mean right bound pressure
Xs = mean(xp);
for t=t1:tf
    pprof = [ tap(2).p(t),tap(3).p(t),tap(4).p(t),tap(5).p(t)];
    Xs = shock_fit(pprof,xp,ps1,ps2,Xs);
    %plot( xp, pprof,'k');hold on;
    %plot( [Xs,Xs], [min(pprof),max(pprof)],'r');
    %hold off;
    %pause;
    
    cc = t-t1+1;
    Xexp(cc) = Xs;
end

% Normalize
xnorm = (Xexp-mean(Xexp)) / 1.78;  % Subtract out the mean and normalize by the throat height (17.8mm (1.78cm))

% Filter out the point to point jumps
jump = abs( xnorm(1:end-1) - xnorm(2:end) );
for i=1:size(xnorm,2)-1
    
    jump = abs( xnorm(i) - xnorm(i+1));
    if ( jump >.75)
        xnorm(i+1) = xnorm(i);
    end
end

dt = 1/srate;
time = linspace(1,size(xnorm,2),size(xnorm,2)) * dt;
time = time * 1000; % Convert to milliseconds

plot( xnorm ) 

ofile = '../../shock_history/exp.mat';
save(ofile,'xnorm','time');


figure(2)
L=10000;
Fs=srate;
hw = hanning(L,'periodic');

for i=2:5

y = tap(i).PSD(1:L);
y = y.*hw;
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
%Y = fftshift(fft(y,NFFT))/L;
Y = (fft(y,NFFT))/L;
f = Fs/2*linspace(0,1,NFFT/2+1);

c(i-1,:) = 2*abs(f'.*Y(1:NFFT/2+1));
% Plot single-sided amplitude spectrum.
%figure(i)
%semilogx(f,2*abs(f'.*Y(1:NFFT/2+1)));hold all; 
%loglog(f,2*abs(f'.*Y(1:NFFT/2+1)));hold all; 
%xlim([f(ceil(size(f,2)/2)) f(end) ])

% Get the variance of p over x
evar(i) = sum((tap(i).p - tap(i).pbar).^2)/max(size(tap(i).p)) / tap(i).pbar^2;


end
%title('Single-Sided Amplitude Spectrum of y(t)')
%xlabel('Frequency (Hz)')
%ylabel('|Y(f)|')


tt = 5.4e-5;
%tt = 3.4e-5;
xoff = 7;

%% Reduce data to a mere N points
N = 1000;
f1 = f(2);
f2 = f(end);
FF = logspace(log10(f1),log10(f2),N);
for i=1:4
    b(i,:) = interp1(f,c(i,:),FF);
    b(i,:) = b(i,:) / sum(b(i,:));
    %b(i,:) = b(i,:) / max(b(i,:));
end

%b = b/N;

% Make a 2d Contour plot of exp. data
lfreq = size(FF,2);
x = zeros(4,lfreq);
y = zeros(4,lfreq);
for i=1:lfreq
    x(:,i) = [0, 1.27, 2.54, 3.81] + 5.5;
    % Filter this term
    %for k=1:4
    %    for ii=2:min(size(b))-1
    %        b(ii,:)=(b(ii-1,:)+b(ii+1,:))/2;
    %    end
    %end
end
Ht = 1.78;
x = x / Ht;

for i=1:4
    y(i,:) = FF*tt; %log10(FF*tt);
    y(i,:) = log10(FF*tt);
    % Filter this term
    for k=1:2
        for ii=2:max(size(b))-1
            b(:,ii)=(b(:,ii-1)+b(:,ii+1))/2;
        end
    end
end
b = b / max(max(b));


%write_visit('../../../matlab_scripts/DATA/exp_spec.tec',x,y,b);

figure(7);
[cp,cp] = contourf(x,y,b,12);
set(cp,'edgecolor','none');
xray = flipud(gray);
colormap(xray)


figure(2);
eX = x(:,1);
plot(eX,evar(2:5),'ko');

figure(2);hold on
load('../sim_var.mat')
sX = X(1:max(size(var)));
plot(sX,var);

save var_comp.mat sX var eX evar;
