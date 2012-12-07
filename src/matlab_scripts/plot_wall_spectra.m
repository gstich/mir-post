clear all;
clc;
close all;


Ht = 1.78;
Up = 32940*1.603;
Fs = 200e3;
LW = 1.5;

% Read in the .tec files and make contour plots
res = 'medInlet';

switch res
    case('corInlet')
        load('DATA/Pw_fullcorInlet.mat');
        dx = 1;
        x1 = 1;
        xn = 470;
        xs = 8;
        %Ptop = Pbot;

    case('ucorInlet')
        load('DATA/Pw_fullucorInlet.mat');
        dx = 1;
        x1 = 1;
        xn = 350;
        xs = 8;

    case('medInlet')
        load('DATA/Pw_fullmedInlet.mat');
        dx = 1;
        x1 = 1;
        xn = 700;
        xs = 8;
end

load(['DATA/sim_',res,'_data.mat'])





%%% Patch on the high wave number decay...
hwn1d = load('../data/Tspec.mat');
kfact = 1/1000*Ht/Up*Fs;
sfact = max(hwn1d.corHWN)*1.5;
f = hwn1d.corHWN(1:200)/sfact;
for i=1:20
    f = gfilter(f);
end
Xhwn = repmat(XX(:,1),1,size(f,1));
kk = kfact * hwn1d.corHWN_k(1:200);
Khwn = repmat(kk,1,size(XX,1))';
Fhwn = repmat(f,1,size(XX,1))';
Fhwn = Fhwn + rand(size(Fhwn))*.15;

contourf(Xhwn,Khwn,Fhwn,32,'edgecolor','none');
hold on;

crange = linspace(min(min(pgrid)),max(max(pgrid))*.5,32);
contourf(XX(:,1:10:end),FF(:,1:10:end),pgrid(:,1:10:end),crange,'edgecolor','none');
xlim([0 6.5]);
ylim([10^-3 5*10^1]);
xray = flipud(gray);
colormap(xray);
set(gca,'YScale','log');



plt.fig=1;
plt.pdf=true;
plt.AR = 2.5;
plt.file = [res,'_2dspectraD'];
plt.xlabel = '$x/H_t$';
plt.ylabel = '$f H_t/U_p$';
pretty_plot(plt);
