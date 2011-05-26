
clc;

% Some parameters for ease of use
X=20; Y=21; U=1; V=2; W=3; RHO=5; P=4; T=6;UU=7;VV=8;WW=9;
MU=16;MUa=17;MUb=18;MUc=19;

Ru = 8.314472e7;
mw = 28.966;
R = Ru/mw;

resolution = 'medium';

switch resolution
    case 'medium',
t1 = 650;
t2 = 675;
Nwall = 300;
sep = 40;
Nst = 8;
for i=1:Nst
    station(i) = Nwall + (i-1)*sep;
end
buff = 5;
path = '/p/lscratchd/olson45/nozzle/post_procmedium3d';
ofile = '../data/BLprof/medium.dat';
    case 'coarse',
t1 = 800;
t2 = 840;
Nwall = 300;
sep = 40;
Nst = 8;
for i=1:Nst
    station(i) = Nwall + (i-1)*sep;
end
buff = 5;
path = '/p/lscratchd/olson45/nozzle/post_proccoarse3d';
ofile = '../data/BLprof/coarse.dat'; 
    case 'fine',
t1 = 195;
t2 = 206;
Nwall = 534;
sep = 50;
Nst = 8;
for i=1:Nst
    station(i) = Nwall + (i-1)*sep;
end
buff = 5;
path = '/p/lscratchd/olson45/nozzle/post_procfine3d';
ofile = '../data/BLprof/fine.dat'; 
side = 'bot';

end

if ~exist('Pmean')
    Pmean = sum_planes(path,'post',t1,t2);
end

if ~exist('fluc')
    fluc = Pmean;
end
    
    
nx = size(Pmean,1);
ny = size(Pmean,2);

%figure(1);hold all;

for ii=1:Nst
    Nwall = station(ii);

% Make a 1D array
Yprof = sum( Pmean(Nwall-buff:Nwall+buff,1:ny/2,:),1);
Yprof2 = sum( Pmean(Nwall-buff:Nwall+buff,ny:-1:ny/2+1,:),1);
Yprof2(:,:,Y) = -Yprof2(:,:,Y);
if(side=='top')
    Yprof = Yprof2;
elseif(side=='bot')
    
else
Yprof = Yprof + Yprof2;
Yprof = Yprof / 2;
end
Yprof = Yprof / (2*buff + 1);


% Thermo-dynamic vars
Pinitial = 1e6;
Tinitial = 300;
NPR = 1.7;
Pin = .5283;
rhoin = .63395;
Mach = .99997;
ReBL = 10e3;
delBL= 2.0e-1;
gamma = 1.4;

%Pin  =   Pmean(1,ny/2, P );
%rhoin  = Pmean(1,ny/2,RHO);
%Uin = Pmean(1,ny/2,U) %Mach*sqrt(Pin*gamma/rhoin);
%mu_0 = Uin * rhoin * delBL / ReBL

rhoin = rhoin*NPR*Pinitial/(Tinitial*R);
Pin = Pin*NPR*Pinitial                    ;
Uin = Mach*sqrt(Pin*gamma/rhoin)         ;
ein = (Pin/(gamma-1))/rhoin            ;
Tin = Pin / (rhoin * R)               ;
mu_0 = Uin * rhoin * delBL / ReBL       
  




% Sutherland's law
T_0 = 273.15;  % Kelvin reference temperature
ST = 110.4;    % Sutherland temperature

T_wall = Yprof(1,1,T)
mu_w = mu_0 * ((T_wall/T_0)^(3.0D0/2.0D0) * (T_0 + ST) / (T_wall + ST) )
mu_w = (Yprof(1,1,MU) + Yprof(1,1,MU) )/ 2    %9.1e-4;

if(ii==1)
% Scaling
rho_w = Yprof(1,1,RHO);
dudy = ( Yprof(1,2,U) - Yprof(1,1,U) )/ ( Yprof(1,2,Y) - Yprof(1,1,Y) );
tauw = mu_w * dudy;
utau = sqrt( tauw / rho_w );
del = mu_w / ( rho_w * utau );
end


figure(1);
y1 = Yprof(1,:,Y)-Yprof(1,1,Y);
u1 = Yprof(1,:,U);

% Van Driest
uvd(1) = 0;
for i=2:ny/2
    dup = u1(i) - u1(i-1);
    uvd(i) = uvd(i-1) + sqrt( Yprof(1,i,RHO) / rho_w) * dup;
end

y1 = y1/del;
up = u1/utau;
uvd = uvd/utau;

xoff = Pmean(station(ii),1,X);
uvd = uvd / max(uvd) + xoff;

%semilogy( uvd,y1,'bo');hold on;

y2 = y1*del/delBL;

figure(1);hold on;
for i=1:ny/2
    x = uvd(i);
    y = y2(i);
    if(x < xoff )
        sym = 'ro';
    else
        sym = 'ko';
    end
    plot( x, y , sym);
end

%x1 = logspace(-1,1.2,20);
%x2 = logspace(.8,3,20);
%k = .41;
%C = 5.1;
hold on;
%semilogx(x1,x1,'k--');hold on;
%semilogx(x2,1/k*log(x2)+C,'k--');
end



