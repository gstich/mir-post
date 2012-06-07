clear all;
clc;


% Time lag computation.
var_map

figure(4);load('../data/fullucorInlet.mat');hold on;
a = sqrt( dataC(:,:,P)*1.4./dataC(:,:,RHO));
M=dataC(:,:,U)./a;
plot(dataC(335,:,Y),M(335,:),'b-');
%disp(trapz(dataC(335,:,Y)-dataC(335,1,Y).*M(335,:))/trapz(dataC(335,:,Y)-dataC(335,1,Y)))
Ubar = sign(dataC(335,:,U));
%Ubar = Ubar.^2;
disp(trapz( Ubar.* (dataC(335,:,Y)-dataC(335,1,Y)))/trapz( max(Ubar).* (dataC(335,:,Y)-dataC(335,1,Y))))


figure(4);load('../data/fullucorAR1.5.mat');hold on;
a = sqrt( dataC(:,:,P)*1.4./dataC(:,:,RHO));
M=dataC(:,:,U)./a;
plot(dataC(335,:,Y),M(335,:),'r');
Ubar = sign(dataC(335,:,U));
%Ubar = Ubar.^2;
disp(trapz( Ubar.* (dataC(335,:,Y)-dataC(335,1,Y)))/trapz( max(Ubar).* (dataC(335,:,Y)-dataC(335,1,Y))))


figure(4);load('../data/fullucorAR1.5NPR1.55.mat');hold on;
a = sqrt( dataC(:,:,P)*1.4./dataC(:,:,RHO));
M=dataC(:,:,U)./a;
plot(dataC(335,:,Y),M(335,:),'r--');
%disp(trapz(dataC(335,:,Y)-dataC(335,1,Y).*M(335,:))/trapz(dataC(335,:,Y)-dataC(335,1,Y)))
Ubar = sign(dataC(335,:,U));
%Ubar = Ubar.^2;
disp(trapz( Ubar.* (dataC(335,:,Y)-dataC(335,1,Y)))/trapz( max(Ubar).* (dataC(335,:,Y)-dataC(335,1,Y))))


figure(4);load('../data/fullucorAR1.7NPR1.70.mat');hold on;
a = sqrt( dataC(:,:,P)*1.4./dataC(:,:,RHO));
M=dataC(:,:,U)./a;
plot(dataC(335,:,Y),M(335,:),'g-');
%disp(trapz(dataC(335,:,Y)-dataC(335,1,Y).*M(335,:))/trapz(dataC(335,:,Y)-dataC(335,1,Y)))
Ubar = sign(dataC(335,:,U));
%Ubar = Ubar.^2;
disp(trapz( Ubar.* (dataC(335,:,Y)-dataC(335,1,Y)))/trapz( max(Ubar).* (dataC(335,:,Y)-dataC(335,1,Y))))


figure(4);load('../data/fullucorAR1.7NPR1.90.mat');hold on;
a = sqrt( dataC(:,:,P)*1.4./dataC(:,:,RHO));
M=dataC(:,:,U)./a;
plot(dataC(335,:,Y),M(335,:),'g--');
%disp(trapz(dataC(335,:,Y)-dataC(335,1,Y).*M(335,:))/trapz(dataC(335,:,Y)-dataC(335,1,Y)))
Ubar = sign(dataC(335,:,U));
%Ubar = Ubar.^2;
disp(trapz( Ubar.* (dataC(335,:,Y)-dataC(335,1,Y)))/trapz( max(Ubar).* (dataC(335,:,Y)-dataC(335,1,Y))))
