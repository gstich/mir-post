% Compare the shock location time history for various cases using the
% pressure signal data.


clc;
close all;
clear all;


Ht = 1.78;
Up = 32940*1.603;
Fs = 200e3;
LW = 1.5;

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

%prob(ii).fname = '../data/shock_history/fullucor.mat';
%prob(ii).side = bot;
%prob(ii).Pcol = 'b';
%prob(ii).range = [1,420,3465-250];
%ii=ii+1;

%prob(ii).fname = '../data/shock_history/fullcor.mat';
%prob(ii).side = top;
%prob(ii).Pcol = 'r';
%prob(ii).range = [1,600,3021-150];
%ii=ii+1;

%prob(ii).fname = '../data/shock_history/fullmed.mat';
%prob(ii).side = top;
%prob(ii).Pcol = 'g';
%prob(ii).range = [1,600,3401-350];
%ii=ii+1;

%prob(ii).fname = '../data/shock_history/fullMU.mat';
%prob(ii).side = top;
%prob(ii).Pcol = 'c';
%prob(ii).range = [1,1,1737];
%ii=ii+1;

%prob(ii).fname = '../data/shock_history/full3xbl.mat';
%prob(ii).side = top;
%prob(ii).Pcol = 'c';
%prob(ii).range = [1,1,1596];
%ii=ii+1;

%prob(ii).fname = '../data/shock_history/fullLam.mat';
%prob(ii).side = top;
%prob(ii).Pcol = 'g';
%prob(ii).range = [1,1,2730];
%ii=ii+1;

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

prob(ii).fname = '../data/shock_history/fullucorInlet.mat';
prob(ii).side = bot;
prob(ii).Pcol = 'b--';
prob(ii).range = [1,1,2600-50];
ii=ii+1;

%prob(ii).fname = '../data/shock_history/fullcorInlet.mat';
%prob(ii).side = top;
%prob(ii).Pcol = 'r--';
%prob(ii).range = [1,1,780-50];
%ii=ii+1;

prob(ii).fname = '../data/shock_history/fullucorAR1.5.mat';
prob(ii).side = bot;
prob(ii).Pcol = 'r-';
prob(ii).range = [1,10,3254];
ii=ii+1;

prob(ii).fname = '../data/shock_history/fullucorAR1.5NPR1.55.mat';
prob(ii).side = bot;
prob(ii).Pcol = 'g-';
prob(ii).range = [1,10,3292];
ii=ii+1;

prob(ii).fname = '../data/shock_history/fullucorAR1.7NPR1.70.mat';
prob(ii).side = bot;
prob(ii).Pcol = 'c-';
prob(ii).range = [1,10,3185];
ii=ii+1;

prob(ii).fname = '../data/shock_history/fullucorAR1.7NPR1.90.mat';
prob(ii).side = bot;
prob(ii).Pcol = 'y-';
prob(ii).range = [1,400,3157];
ii=ii+1;



expS(1).fname = '../data/shock_history/exp.mat';
expS(1).side = top;
expS(1).Pcol = 'k--';
xFs = 200e3;

% PDF x sample space
xx = linspace(-2,2,400);
ssig = .1^2;

for i=1:size(prob,2)
	disp(prob(i).fname)
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




for i=1:0 %size(expS,2)

    nbin = 20;
    over = .5;
    
    load(expS(i).fname);
    nT = max(size(xnorm));
    lim = nT/nbin;
    %Xs = xnorm';
    Xs = xnorm(1:lim)';
    timeT = time(1:lim);
    mean(Xs)
    Xs = (Xs - mean(Xs))/Ht;
    FXs = compute_spectra(Xs);
    FXs(1,:) = FXs(1,:) * xFs *Ht/Up;
    
    
    FXs(2,:) = 0;

    % For loop to average bins

    for j=1:nbin/over - 1
      Xs = xnorm( (j*over-over)*lim+1:(j*over-over)*lim + lim  )';
      Xs = (Xs-mean(Xs))/Ht;
      FX = compute_spectra(Xs);
      FXs(2,:) = FXs(2,:) + FX(2,:);
    end
    %FXs(1,:) = FXs(1,:)*nbin;
    FXs(2,:) = FXs(2,:) / (nbin/over - 1);  % Average over bins
    FXs(2,:) = FXs(2,:) * (xFs * (lim/xFs)^.5  ) *Ht/Up;



    
    
    % Spectra
    figure(1)
    loglog(FXs(1,:),FXs(2,:),expS(i).Pcol,'Linewidth',2);
    hold all
    
    % Compensated spectra
    figure(2)
    semilogx(FXs(1,:),FXs(2,:).*FXs(1,:),'k*');
    hold all
    
    % Raw data  (just of last bin)
    figure(3)
    plot(timeT,Xs,expS(i).Pcol);
    hold all
    
    % PDF (rms f(x))
    %Xs = xnorm;
    %Xs = (Xs - mean(Xs))/Ht;
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


figure(1)
