clear all;
clc;
close all;

%% Plot the BL data: Mean and Reynolds stress
pdfE = true;
%pdfE = false;
rez = 'fine';

Ht = 1.78;

res(1).name = rez;
ref(1) = 1;

% Figure name
figs(1).name = 'mean_st';

% Plot symbols
res(1).sym = 'k-';

% Figure option
LW = 2;         % LineWidth
FSn = 25;       % FontSize labels
FSa = 18;       % FontSize axis

% Data format key
key1 = '%% < 1-5 > y(cm), y+, del, u(cm/s),u+';
key2 = '%% < 6-11> uvd, utau, uu, vv, ww, delBL';
key.y = 1;      key.yp = 2;     key.del = 3;
key.u = 4;      key.up = 5;     key.uvd = 6;
key.utau = 7;   key.x  = 8;     key.delBL = 9; 
key.tauw = 10;

% Load the data
for i=1:size(res,2)
    file = ['BLprof/',res(i).name, '_st.mat'];
    load(file);
    res(i).data = pprof;
end


x1 = logspace(-1,1.2,20);
x2 = logspace(.8,3,20);
k = .37;
C = 5.2;
%k = .41;
%C = 5.1;


ss = 1;
for ss=1:size(pprof,3)
    
for i=1:size(res,2)
    
    % Get reference utau,del,delBL,tau_w
    utau = res(ref(i)).data(1,key.utau,ss);
    del = res(ref(i)).data(1,key.del,ss);
    delBL = res(ref(i)).data(1,key.delBL,ss);
    tauw = res(ref(i)).data(1,key.tauw,ss);
    
    y = res(i).data(:,key.y,ss) / del; y = y - y(1);
    u = res(i).data(:,key.u,ss) / utau;
    utau_i = res(i).data(1,key.utau,ss);
    uvd = res(i).data(:,key.uvd,ss);
    uvd = uvd * utau_i / utau;
    semilogx(y,uvd,'LineWidth',LW);hold all;
end
end


semilogx(x1,x1,'k--');hold on;
semilogx(x2,1/k*log(x2)+C,'k--');hold on;



h1 = xlabel(['$y^+$']);
set(h1,'Interpreter','latex','FontSize',FSn);
h2 = ylabel('$U^+_{VD}$');
set(h2,'Interpreter','latex','FontSize',FSn);
box on;
set(gca,'FontSize',FSa);
xlim([.5 .5e3]);
ylim([0 25]);
set(gca,'Position',[.13,.13,.775,.8107])



%% Add equations to figure(1)
figure(1);
% Add the text
tx = 2.2;
ty = 12;
h1 = text(tx,ty,'$U^+_{VD} = y^+$');
h2 = text(tx+25,ty,'$U^+_{VD} = \log(y^+)/.37+5.2$');
set(h1, 'interpreter', 'latex','FontSize',FSa)
set(h2, 'interpreter', 'latex','FontSize',FSa,'rotation',30)

for i=1:4
L(i).n = ['$x/H_t=',num2str(pprof(1,key.x,i)/Ht,2),' $'];
end

h = legend(L(1).n,L(2).n,L(3).n,L(4).n);
set(h, 'interpreter', 'latex','FontSize',FSa,'Location','NorthWest');
legend boxoff;


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

