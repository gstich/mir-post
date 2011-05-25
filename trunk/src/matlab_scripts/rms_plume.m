clear all;
clc;
close all;


path = 'static/';
base = 'pressure_0';

path = '../visit_scripts/';
base = 'ptot_';

ext = '.dat';
xnum = 1;
ynum = 2;
znum = 3;
Vnum = 4;

npts = 201;
t0 = 560;
tf = 595;
p0 = 1e6;
NPR = 1.7;
pres = p0*NPR;
Ht = 1.78;

Vdata = zeros(npts,tf-t0+1,2);

% Read the entire data-set by looping through the files
for i=t0:tf
    file = [ path, base, num2str(i,'%4.4i'),ext];
    disp(file);
    d = load(file);
    
    Vdata(:,i-t0+1,1) = d(:,ynum);
    Vdata(:,i-t0+1,2) = d(:,Vnum);
end


% Get the average pressure at each point
pave = zeros(npts,1);
for i=1:tf-t0+1
    tmp = Vdata(:,i,2);
    pave = pave + tmp; 
end
pave = pave / (tf-t0+1);


% Get the RMS pressure
prms = zeros(npts,1);
for i=1:tf-t0+1
    tmp = Vdata(:,i,2);
    prms = prms + (tmp-pave).^2; 
end
prms = prms / (tf-t0+1);
prms = sqrt(prms);

y = Vdata(:,1,1);
y = y/Ht;
figure(1);
plot(y,prms/pres);



% Read in the exp. data file

% Data extents on plot
xx = [-2,2];
yy = [0.0,0.6];

% Data extents on picture-in pixels
LL = [103,471];
UR = [668, 20];

figure(3);
exp = imread('ptot_rms2.png');

%imshow(exp);
image(exp);
xpix = size(exp,2);
ypix = size(exp,1);


xdata = y;
ydata = prms/pres;

xdata = (xdata - xx(1) ) / (xx(2)-xx(1));
xdata = xdata * (UR(1)-LL(1)) + LL(1);

ydata = -(ydata - yy(2)) / (yy(2)-yy(1));
ydata = ydata * (LL(2)-UR(2)) + UR(2);

hold on;
plot(xdata,ydata,'b--','LineWidth',2);



