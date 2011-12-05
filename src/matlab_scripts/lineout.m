function [v,x,y] = lineout(xx,yy,var,xpt,ypt,Npt)

x = linspace(xpt(1),xpt(2),Npt);
y = linspace(ypt(1),ypt(2),Npt);

v = griddata(xx,yy,var,x,y);


end