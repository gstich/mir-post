clear all;
clc;
close all;

%% Plot the BL data: Mean and Reynolds stress
%pdfE = true;
pdfE = false;
res(1).name = 'coarse';
res(2).name = 'medium';
res(3).name = 'fine';
ref(1) = 1;
ref(2) = 2;
ref(3) = 2;

% Figure name
figs(1).name = 'mean';
figs(2).name = 'rms';

% Plot symbols
res(1).sym = 'k--';
res(2).sym = 'k-.';
res(3).sym = 'k-';

% Figure option
LW = 2;         % LineWidth
FSn = 25;       % FontSize labels
FSa = 18;       % FontSize axis

% Data format key
key1 = '%% < 1-5 > y(cm), y+, del, u(cm/s),u+';
key2 = '%% < 6-11> uvd, utau, uu, vv, ww, delBL';
key.y = 1;      key.yp = 2;     key.del = 3;
key.u = 4;      key.up = 5;     key.uvd = 6;
key.utau = 7;   key.uu = 8;     key.vv = 9;
key.ww = 10;    key.delBL = 11; key.tauw = 12;

% Load the data
for i=1:size(res,2)
    file = ['BLprof/',res(i).name, '.dat'];
    res(i).data = load(file);
end


x1 = logspace(-1,1.2,20);
x2 = logspace(.8,3,20);
%k = .375;
%C = 6.4;
k = .41;
C = 5.1;

%% Mean Profiles
figure(1);clf;
semilogx(x1,x1,'k--');hold on;
semilogx(x2,1/k*log(x2)+C,'k--');hold on;


for i=1:size(res,2)
    
    % Get reference utau,del,delBL,tau_w
    utau = res(ref(i)).data(1,key.utau);
    del = res(ref(i)).data(1,key.del);
    delBL = res(ref(i)).data(1,key.delBL);
    tauw = res(ref(i)).data(1,key.tauw);
    
    y = res(i).data(:,key.y) / del; y = y - y(1);
    u = res(i).data(:,key.u) / utau;
    utau_i = res(i).data(1,key.utau);
    uvd = res(i).data(:,key.uvd);
    uvd = uvd * utau_i / utau;
    semilogx(y,uvd,res(i).sym,'LineWidth',LW);hold on;
end
figure(1);


h1 = xlabel(['$y^+$']);
set(h1,'Interpreter','latex','FontSize',FSn);
h2 = ylabel('$u^+$');
set(h2,'Interpreter','latex','FontSize',FSn);
box on;
set(gca,'FontSize',FSa);
xlim([1 2e3]);
ylim([0 25]);
set(gca,'Position',[.13,.13,.775,.8107])

figure(2);clf;

for i=1:size(res,2)
    y = res(i).data(:,key.y) / delBL; y = y - y(1);
    uu = sqrt( res(i).data(:,key.uu) / utau^2 );    
    vv = sqrt( res(i).data(:,key.vv) / utau^2 );
    ww = sqrt( res(i).data(:,key.ww) / utau^2 );
    plot(y,uu,res(i).sym,y,vv,res(i).sym,y,ww,res(i).sym,'LineWidth',LW);hold on;
end

figure(2);

h1 = xlabel(['$y/\delta$']);
set(h1,'Interpreter','latex','FontSize',FSn);
h2 = ylabel('$\sqrt{\overline{u^\prime u^\prime}_i}/u_\tau$');
set(h2,'Interpreter','latex','FontSize',FSn);
set(gca,'FontSize',FSa);
box on;
set(gca,'FontSize',FSa);
xlim([0 1.2])
ylim([0 3.5])
set(gca,'Position',[.13,.13,.775,.8107])


% Load and plot Spalart DNS
%dns = load('BLprof/spalart_dns.dat');
%figure(1);hold on;
%plot( 10.^dns(:,1), dns(:,2),'k-');


% Save the figures and convert them to .pdf
if (pdfE)
    for i=1 : size (figs , 2)
        fname = [ '../figs/',figs(i).name , '.eps' ];
        figure(i);
        print('-depsc2',fname)
        eps2pdf(fname)
        delete(fname)
    end
end

