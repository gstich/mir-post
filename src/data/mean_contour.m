
function mean_contour(x,y,var,xrng,yrng,vrng,ncon,iFig,xL,yL,cbL)

% Figure option
LW = 2;         % LineWidth
FSn = 25;       % FontSize labels
FSa = 18;       % FontSize axis

AR = (xrng(2)-xrng(1)) / (yrng(2)-yrng(1));

figure(iFig);

if(vrng(1)==vrng(1))
    Ncon = ncon;
    vrng(1) = min(min(var));
    vrng(2) = max(max(var));
else
    Ncon = linspace(vrng(1),vrng(2),ncon);
end

contourf(x,y,min(max(var,vrng(1)),vrng(2)),Ncon,'EdgeColor','none');
axis equal; xlim(xrng); ylim(yrng); hold on; 
cb = colorbar; h2 = ylabel(cb,cbL);set(h2,'Interpreter','latex','FontSize',FSn);
h = xlabel(xL); set(h,'Interpreter','latex','FontSize',FSn);
h1 = ylabel(yL); set(h1,'Interpreter','latex','FontSize',FSn);
set(gca,'FontSize',FSa);


g = figure(iFig);
set(g,'Position',[200,395,250*AR+50,250+50])
set(g,'PaperSize',[8,8/AR]);
set(g,'PaperPositionMode','auto');

%% Odd and even for prob location
xoff = 5;
if(mod(iFig,2)==1)
    ix = [380,87] - xoff;
    iy = [80,18];
else
    ix = [577,131] - xoff;
    iy = [95,21];
end

hold on;
plot(x(ix(1),iy(1)),y(ix(1),iy(1)),'go','LineWidth',2,'MarkerFaceColor','k','MarkerSize',15);hold on
plot(x(ix(2),iy(2)),y(ix(2),iy(2)),'go','LineWidth',2,'MarkerFaceColor','k','MarkerSize',15);





end