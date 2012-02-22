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

% Figure option
LW = 1;         % LineWidth
FSn = 25;       % FontSize labels
FSa = 18;       % FontSize axis

pdfE = false;
pdfE = true;

f(1).name = 'XSspectraT';
f(2).name = 'XScompT';
f(3).name = 'XSrawT';
f(4).name = 'XSvarT';

%path = '~/Dropbox/Britton/THESIS/Figures/nozzle/convergence/';
path = '/Users/bolson/Dropbox/Britton/THESIS/Figures/nozzle/convergence/';


ii=1;
% prob(ii).fname = '../data/shock_history/coarse.mat';
% prob(ii).side = top;
% prob(ii).Pcol = 'k';
% ii=ii+1;

% prob(ii).fname = '../data/shock_history/medium.mat';
% prob(ii).side = bot;
% prob(ii).Pcol = 'b';
% ii=ii+1;


prob(ii).fname = '../data/shock_history/coarse.mat';
prob(ii).side = top;
prob(ii).Pcol = 'k';
%prob(ii).range = [1,1,4445];
prob(ii).range = [1,800,2946];
ii=ii+1;


prob(ii).fname = '../data/shock_history/mediumv2.mat';
prob(ii).side = mid;
prob(ii).Pcol = 'b';
%prob(ii).range = [1,500,1500];
prob(ii).range = [1,11,2254];
ii=ii+1;

%prob(ii).fname = '../data/shock_history/mediumNarrow.mat';
%prob(ii).side = mid;
%prob(ii).Pcol = 'c';
%prob(ii).range = [1,1,1123];
%ii = ii + 1;
 
%prob(ii).fname = '../data/shock_history/fullwide.mat';
%prob(ii).side = mid;
%prob(ii).Pcol = 'g';
%prob(ii).range = [1,200-150,1151-150];
%ii = ii + 1;




expS(1).fname = '../data/shock_history/exp.mat';
expS(1).side = top;
expS(1).Pcol = 'r--';
xFs = 200e3;

% PDF x sample space
xx = linspace(-2,2,400);
ssig = .1^2;

for i=1: 1 %size(prob,2)

    load(prob(i).fname);
    Xs = XSS(:,prob(i).side);
    if (prob(i).range(1) == 1)
        i1 = prob(i).range(2);
        iN = prob(i).range(3);
        Xs = XSS(i1:iN,prob(i).side);
        tt = time;
        clear time
        time = tt(i1:iN);
        time = time - time(1);
    end
    mean(Xs)
    Xs = (Xs - mean(Xs))/Ht;
    FXs = compute_spectra(Xs);
    FXs(1,:) = FXs(1,:) * Fs *Ht/Up;
    
    
    % Spectra
    figure(1)
    loglog(FXs(1,:),FXs(2,:),prob(i).Pcol,'Linewidth',LW);
    hold all
    
    % Compensated spectra
    figure(2)
    semilogx(FXs(1,:),FXs(2,:).*FXs(1,:)*Up/Ht,prob(i).Pcol,'Linewidth',LW);
    hold all
    
    % Raw data
    figure(3)
    plot(time*Up/Ht/1000,Xs,prob(i).Pcol,'Linewidth',LW);
    hold all
    
    
    % PDF (rms f(x))
    w = xx * 0.0;
    for j = 1:size(Xs,1)
        w = w + exp( -(xx-Xs(j)).^2/ssig );
    end
    w = w/size(Xs,1);
    
    figure(4);
    plot(xx,w,prob(i).Pcol,'Linewidth',LW*2);
    hold all;
    
end




for i=1:size(expS,2)

    lim = 20000;
    load(expS(i).fname);
    %Xs = xnorm';
    Xs = xnorm(1:lim)';
    timeT = time(1:lim);
    Xs = (Xs - mean(Xs))/Ht;
    FXs = compute_spectra(Xs);
    FXs(1,:) = FXs(1,:) * xFs *Ht/Up;
    
    
    % Spectra
    figure(1)
    loglog(FXs(1,:),FXs(2,:),expS(i).Pcol,'Linewidth',LW);
    hold all
    
    % Compensated spectra
    figure(2)
    semilogx(FXs(1,:),FXs(2,:).*FXs(1,:)*Up/Ht,expS(i).Pcol,'Linewidth',LW);
    hold all
    
    % Raw data
    figure(3)
    plot(timeT*Up/Ht/1000,Xs,expS(i).Pcol,'Linewidth',LW);
    hold all
    
    % PDF (rms f(x))
    w = xx * 0.0;
    for j = 1:size(Xs,1)
        w = w + exp( -(xx-Xs(j)).^2/ssig );
    end
    w = w/size(Xs,1);
    
    figure(4);
    plot(xx,w,expS(i).Pcol,'Linewidth',LW*2);
    hold all;
    
end


figure(3);
xlim([0 20])
