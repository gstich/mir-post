clear all;
clc;
close all;


% Compare mean solutions to establish grid convergence.


% Load the 2 mesh sizes
load('../data/2dmv2.mat');
load('../data/2dcv2.mat');
Ht = 1.78;
Up = 32940*1.603;

% Get some variables
var_map
uc = dataC(:,:,U);
um = dataM(:,:,U);
xc = dataC(:,:,X);yc = dataC(:,:,Y);
xm = dataM(:,:,X);ym = dataM(:,:,Y);


% Contour plot
figure(1);
cnts = linspace(-1e4,6e4,16);
contour(xc,yc,uc,cnts,'color','black');
hold on;
contour(xm,-ym,um,cnts,'color','blue');
legend('Coarse','Medium')
xlim([-2 17]);
ylim([-2 2]);


% Line out of data.
clear urkm
clear urkc
xLO = 17;
epsLO = .2;


% % Coarse line out
% ii=1;
% for i=ceil(size(xc,1)/2):size(xc,1)
%     for j=1:size(xc,2)
%         if(xc(i,j) > xLO-epsLO && xc(i,j) < xLO+epsLO)
%             urkc(1,ii)=yc(i,j);
%             urkc(2,ii)=uc(i,j);
%             ii=ii+1;
%         end
%     end
% end
% [urkc(1,:),I] = sort(urkc(1,:));
% urkc(2,:) = urkc(2,I); 
% 
% % Medium line out
% ii=1;
% for i=ceil(size(xm,1)/2):size(xm,1);
%     for j=1:size(xm,2);
%         if(xm(i,j) > xLO-epsLO && xm(i,j) < xLO+epsLO);
%             urkm(1,ii)=ym(i,j);
%             urkm(2,ii)=um(i,j);
%             ii=ii+1;
%         end;
%     end;
% end
% [urkm(1,:),I] = sort(urkm(1,:));
% urkm(2,:) = urkm(2,I);

% Interpolate onto line med
NN = 25;
[ZIc,XI,YIc] = lineout(xc(ceil(end/2):end,:),yc(ceil(end/2):end,:),uc(ceil(end/2):end,:),[xLO,xLO],[-3,3]*Ht,NN);

[ZIm,XI,YIm]= lineout(xm(ceil(end/2):end,:),ym(ceil(end/2):end,:),um(ceil(end/2):end,:),[xLO,xLO],[-3,3]*Ht,NN);


figure(2);
plot(YIc/Ht,ZIc/Up,'k-','LineWidth',3);
hold on;
plot(-YIm/Ht,ZIm/Up,'b-','LineWidth',3);
legend('Coarse','Medium')

xlim([-3 3]);
