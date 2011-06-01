function shock_track

clear all;
clc;
close all;

%% Load the pressure data on interval... average and plot



pref = 7.5e5;
dt = 5.0e-6;


data = '/p/lscratchd/olson45/nozzle/';
dir = [data,'nozzlecoarse3d/'];
pa = zeros(512,4);
intv = [1500,1916];

%dir = [data,'nozzlemedium3d/'];
%pa = zeros(768,4);
%intv = [800,1060];

%weight = abs(intv(1)-intv(2)) + 1;

ii = 0;
for i=intv(1):intv(2)
    ii = ii + 1;
    file = [ dir ,'vis',num2str(i,'%4.4i'),'/pressure.dat' ]
    
    pf = load(file);
    
    pa = pa + pf;
    
    x = pf(:,1);
    bot = pf(:,2);
    top = pf(:,3);
    mid = pf(:,4);
    
    val = top;
    [xs,ps] = get_shock_location(x,val);
    XSS(ii,1) = xs;
    
    %val = mid;
    %[xs,ps] = get_shock_location(x,val);
    XSS(ii,2) = xs;
    
    %val = bot;
    %[xs,ps] = get_shock_location(x,val);
    XSS(ii,3) = xs;
    
    
    
    %figure(3);
    %plot(x,val);hold on;
    %plot(X1,Y1,'go');%plot(X2,Y2,'go');plot(X3,Y3,'go');
    %plot(xs(ii),ps(ii),'ro');drawnow;pause(.01);hold off;
    %hold off;
    
end

%Xs(:,1) = xs;
%for j=1:5
%Xs = gfilter(Xs);
%end
%plot(Xs);

steps = intv(2)-intv(1)+1;
time = linspace(0,(steps-1)*dt,steps);
time = time * 1000;  % make mili-seconds


plot(time,XSS(:,1),time,XSS(:,2),time,XSS(:,3));
%xmean = mean(Xs);




% Read in the exp. data file

% Data extents on plot
xx = [0,4];
yy = [-0.5,0.7];

% Data extents on picture-in pixels
LL = [509,348];
UR = [712, 14];

%figure(3);
%exp = imread('shock_hist.png');

%imshow(exp);
% xpix = size(exp,2);
% ypix = size(exp,1);
% 
% 
% xdata = time;
% ydata = (Xs-xmean)/2.23;
% 
% xdata = (xdata - xx(1) ) / (xx(2)-xx(1));
% xdata = xdata * (UR(1)-LL(1)) + LL(1);
% 
% ydata = -(ydata - yy(2)) / (yy(2)-yy(1));
% ydata = ydata * (LL(2)-UR(2)) + UR(2);
% 
% hold on;
% plot(xdata,ydata,'k-','LineWidth',1.5);

xs = XSS(:,1);
save Xs;




end

function [xs,ps] = get_shock_location(x,press)

    val = press;
    for j=1:50
        val = gfilter(val);
    end
    
    
    %% Try fitting parabola to min.
    off = 1;
    [Pmin,Nmin] = min(val);
    X1 = x(Nmin);
    Y1 = Pmin;
    X2 = x(Nmin - off);
    Y2 = val(Nmin - off);
    X3 = x(Nmin + off);
    Y3 = val(Nmin + off);
    
    DD = X1*(X3^2-X2^2)-X2*X3^2 + X2^2*X3 + X1^2*(X2-X3);
    A = X1*(Y3-Y2)-X2*Y3+X3*Y2+(X2-X3)*Y1;
    A = A/DD;
    B = X1^2*(Y3-Y2)-X2^2*Y3+X3^2*Y2+(X2^2-X3^2)*Y1;
    B = B/DD;
    
    xs =  B / (2*A);
    ps = Pmin;

end



