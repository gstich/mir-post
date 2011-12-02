% Compare the shock location time history for various cases using the
% pressure signal data.


clc;
close all;
clear all;


Ht = 1.78;
Up = 32940*1.603;
Fs = 200e3;

top = 1;
mid = 2;
bot = 3;


ii=1;
% prob(ii).fname = '../data/shock_history/coarse.mat';
% prob(ii).side = top;
% prob(ii).Pcol = 'k';
% ii=ii+1;

% prob(ii).fname = '../data/shock_history/medium.mat';
% prob(ii).side = bot;
% prob(ii).Pcol = 'b';
% ii=ii+1;


prob(ii).fname = '../data/shock_history/coarsev2.mat';
prob(ii).side = top;
prob(ii).Pcol = 'r';
ii=ii+1;


prob(ii).fname = '../data/shock_history/mediumv2.mat';
prob(ii).side = mid;
prob(ii).Pcol = 'c';

expS(1).fname = '../data/shock_history/exp.mat';
expS(1).side = top;
expS(1).Pcol = 'g';
xFs = 200e3;

% PDF x sample space
xx = linspace(-2,2,400);
ssig = .1^2;

for i=1:size(prob,2)

    load(prob(i).fname);
    Xs = XSS(:,prob(i).side);
    %if ( i == 1)
    %    Xs = XSS(end/2:end,prob(i).side);
    %    tt = time;
    %    clear time
    %    time = tt(end/2:end);
    %end
    mean(Xs)
    Xs = (Xs - mean(Xs))/Ht;
    FXs = compute_spectra(Xs);
    FXs(1,:) = FXs(1,:) * Fs *Ht/Up;
    
    
    % Spectra
    figure(1)
    loglog(FXs(1,:),FXs(2,:),prob(i).Pcol);
    hold all
    
    % Compensated spectra
    figure(2)
    semilogx(FXs(1,:),FXs(2,:).*FXs(1,:),prob(i).Pcol);
    hold all
    
    % Raw data
    figure(3)
    plot(time,Xs,prob(i).Pcol);
    hold all
    
    
    % PDF (rms f(x))
    w = xx * 0.0;
    for j = 1:size(Xs,1)
        w = w + exp( -(xx-Xs(j)).^2/ssig );
    end
    w = w/size(Xs,1);
    
    figure(4);
    plot(xx,w,prob(i).Pcol);
    hold all;
    
end




for i=1:size(expS,2)

    lim = 3000;
    load(expS(i).fname);
    %Xs = xnorm';
    Xs = xnorm(1:lim)';
    timeT = time(1:lim);
    Xs = (Xs - mean(Xs))/Ht;
    FXs = compute_spectra(Xs);
    FXs(1,:) = FXs(1,:) * xFs *Ht/Up;
    
    
    % Spectra
    figure(1)
    loglog(FXs(1,:),FXs(2,:),expS(i).Pcol);
    hold all
    
    % Compensated spectra
    figure(2)
    semilogx(FXs(1,:),FXs(2,:).*FXs(1,:),expS(i).Pcol);
    hold all
    
    % Raw data
    figure(3)
    plot(timeT,Xs,expS(i).Pcol);
    hold all
    
    % PDF (rms f(x))
    w = xx * 0.0;
    for j = 1:size(Xs,1)
        w = w + exp( -(xx-Xs(j)).^2/ssig );
    end
    w = w/size(Xs,1);
    
    figure(4);
    plot(xx,w,expS(i).Pcol);
    hold all;
    
end


figure(3);
xlim([0 20])