clear all;
clc;


% Time lag computation.
var_map
LW = 1.2;
Ht = 1.78;

FSn = 25;
FSa = 20;

data(1).name = '../data/fullucorInlet.mat'; data(1).col = 'b';
data(2).name = '../data/fullucorAR1.5.mat'; data(2).col = 'r--';
data(3).name = '../data/fullucorAR1.5NPR1.55.mat'; data(3).col = 'r-';
data(4).name = '../data/fullucorAR1.7NPR1.70.mat'; data(4).col = 'g--';
data(5).name = '../data/fullucorAR1.7NPR1.90.mat'; data(5).col = 'g-';

for i=1:5
    
    figure(6);load(data(i).name);hold on;
    a = sqrt( dataC(:,:,P)*1.4./dataC(:,:,RHO));
    M=dataC(:,:,U)./a;
    plot(dataC(335,:,Y)/2/max(dataC(335,:,Y)),M(335,:)/max(M(335,:)) ,data(i).col,'Linewidth',LW);
    
    %disp(trapz(dataC(335,:,Y)-dataC(335,1,Y).*M(335,:))/trapz(dataC(335,:,Y)-dataC(335,1,Y)))
    Ubar = sign(dataC(335,:,U));
    %Ubar = Ubar.^2;
    disp(trapz( Ubar.* (dataC(335,:,Y)-dataC(335,1,Y)))/trapz( max(Ubar).* (dataC(335,:,Y)-dataC(335,1,Y))))
    
    
    % Contour plot of separation region
    figure(i);
    contourf( dataC(:,:,Y)/Ht, dataC(:,:,X)/Ht, M ,32, 'edgecolor','none'); hold on
    contour( dataC(:,:,Y)/Ht, dataC(:,:,X)/Ht, M ,[0 0], 'g--','Linewidth',3.5);
    
    set(gca,'YDir','reverse');
    axis equal;
    xlim([-1 1]);
    ylim([2 12]/Ht);
    
    
    % Set the xlabel
    h = xlabel( '$y/H_t$');
    set(h,'Interpreter','latex','FontSize',FSn);
    
    % Set the ylabel
    h = ylabel( '$x/H_t$');
    set(h,'Interpreter','latex','FontSize',FSn);
    
    
    % Set the axis to the right size
    set(gca,'FontSize',FSa);
    
    box on;
    
    plt1.fig = i;
    %plt1.pdf = false;
    plt1.AR = .5;
    plt1.wide = 300;
    plt1.file = ['mean_Mach_2d_case',int2str(i)];
    plt1.zbuf   = 0.25;    % Buffer on left/bottom of plot between axis/edge
    plt1.fbuf  = .05;
    %pretty_plot(plt1);
    
    
    % Plot the mean pressure
    figure(7);hold on;
    plot( dataC(:,end/2,X), dataC(:,end,Y)/ dataC(1,end,Y),data(i).col);
    
    
end


figure(6);hold on;plot([-.5,.5],[0,0],'k--','Linewidth',1.2);
lkey = {'Case 1','Case 2','Case 3','Case 4','Case 5','Location','NorthWest'};
plt.fig = 6;
%plt.pdf = false;
plt.file = 'mean_U_exit';
plt.xlabel = '$y/H_{\textrm{exit}}$';
plt.ylabel = '$U/U_{\textrm{exit}}$';
plt.legend = lkey;
plt.misc = 'legend boxoff; box on;';
pretty_plot(plt);


% plt1.fig = 1;
% plt1.pdf = false;
% plt1.AR = 4;
% plt1.wide = 1500;
% plt1.file = 'mean_Mach_2d';
% if plt1.pdf
%     fname = [ plt1.file , '.eps' ];
%     print('-depsc2',fname)
%     eps2pdf(fname)
%     delete(fname)
% end

