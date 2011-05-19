clear all;
clc;
close all;

pref = 7.5e5;
dt = 5.0e-6;

%% Load the data
load shock_history/Xs_cor;
Xs_cor = XSS;

load shock_history/Xs_med;
Xs_med = XSS;

steps = size(Xs_cor,1);
timeC = linspace(0,(steps-1)*dt,steps);
timeC = timeC * 1000;

steps = size(Xs_med,1);
timeM = linspace(0,(steps-1)*dt,steps);
timeM = timeM * 1000;

% Figure option
LW = 2;         % LineWidth
FSn = 18;       % FontSize labels
FSa = 12;       % FontSize axis

s = 200;
f = size(Xs_cor,1);

figure(1);hold on;
plot(timeC(1:f-s+1),Xs_cor(s:f,1),timeC(1:f-s+1),Xs_cor(s:f,2),timeC(1:f-s+1),Xs_cor(s:f,3));
plot(timeM,Xs_med(:,1),timeM,Xs_med(:,2),timeM,Xs_med(:,3));
