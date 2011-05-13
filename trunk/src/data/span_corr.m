clear all;
clc;
close all;


%% Plot the BL data: Mean and Reynolds stress
pdfE = true;
res(1).name = 'coarse';
res(2).name = 'medium';
ref = 2;


% Plot symbols
res(1).sym = '--';
res(2).sym = '-.';

% Figure option
LW = 2;         % LineWidth
FSn = 18;       % FontSize labels
FSa = 12;       % FontSize axis
cmax = .2;      % Contour Level Max
Nc = 64;        % Number of contour level

%% Coarse
figure(1);
load('span_corr/cor.mat');
contourf(XX,YY+2,tay,[0:.01:.2],'edgecolor','none')


%% Correlation plots
figure(2);hold on;
x = [1:1:size(Cor,1)];

plot(  x/x(end) *.5, Cor(:,1) ,'--');
plot(  x/x(end) *.5, Cor(:,106) ,'--');
plot(  x/x(end) *.5, Cor(:,399), '--');

%% Medium
figure(1);hold on;
load('span_corr/med.mat');
contourf(XX,YY-2,tay,[0:cmax/Nc:cmax],'edgecolor','none')
colorbar;
axis equal

%% Correlation plots
figure(2);hold on;
x = [1:1:size(Cor,1)];
plot(  x/x(end) *.5, Cor(:,1) );
plot(  x/x(end) *.5, Cor(:,106) );
plot(  x/x(end) *.5, Cor(:,237) );





