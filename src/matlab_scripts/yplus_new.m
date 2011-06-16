
clc;

% Some parameters for ease of use
X=20; Y=21; U=1; V=2; W=3; RHO=5; P=4; T=6;UU=7;VV=8;WW=9;
MU=16;MUa=17;MUb=18;MUc=19;

Ru = 8.314472e7;
mw = 28.966;
R = Ru/mw;

resolution = 'fine';

switch resolution
    case 'medium',
t1 = 700;  %550
t2 = 1100;  %580
Nwall = 135;
buff = 5;
path = '/p/lscratchd/olson45/nozzle/post_procmedium3d';
ofile = '../data/BLprof/medium.dat';
    case 'coarse',
t1 = 1000;
t2 = 1400;
Nwall = 87;
buff = 5;
path = '/p/lscratchd/olson45/nozzle/post_proccoarse3d';
ofile = '../data/BLprof/coarse.dat'; 
    case 'fine',
t1 = 195;
t2 = 233; %205
Nwall = 174;
buff = 5;
path = '/p/lscratchd/olson45/nozzle/post_procfine3d';
ofile = '../data/BLprof/fine.dat'; 

end

if ~exist('Pmean')
    Pmean = sum_planes(path,'post',t1,t2);
end

if ~exist('fluc')
    fluc = Pmean;
end
    
    
nx = size(Pmean,1);
ny = size(Pmean,2);


% Make a 1D array
Yprof = sum( Pmean(Nwall-buff:Nwall+buff,1:ny/2,:),1);
Yprof2 = sum( Pmean(Nwall-buff:Nwall+buff,ny:-1:ny/2+1,:),1);
Yprof2(:,:,Y) = -Yprof2(:,:,Y);
Yprof = Yprof + Yprof2;
Yprof = Yprof / 2;
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

% Scaling
rho_w = Yprof(1,1,RHO);
dudy = ( Yprof(1,2,U) - Yprof(1,1,U) )/ ( Yprof(1,2,Y) - Yprof(1,1,Y) );
tauw = mu_w * dudy;
utau = sqrt( tauw / rho_w );
del = mu_w / ( rho_w * utau );

figure(1);clf;
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

semilogx( y1,uvd,'bo');
x1 = logspace(-1,1.2,20);
x2 = logspace(.8,3,20);
k = .41;
C = 5.1;
hold on;
semilogx(x1,x1,'k--');hold on;
semilogx(x2,1/k*log(x2)+C,'k--');


% Make a 1D array
Fprof = sum( fluc(Nwall-buff:Nwall+buff,1:ny/2,:),1);
Fprof2 = sum( fluc(Nwall-buff:Nwall+buff,ny:-1:ny/2+1,:),1);
Fprof2(:,:,Y) = -Fprof2(:,:,Y);
Fprof = Fprof + Fprof2;
Fprof = Fprof / 2;
Fprof = Fprof / (2*buff + 1);



figure(2);hold on;
y2 = y1*del/delBL;
u1 = Fprof(1,:,UU).*Yprof(1,:,RHO)/tauw;
u1 = sqrt(u1);
plot(y2,u1,'b');hold on;
u1 = Fprof(1,:,VV).*Yprof(1,:,RHO)/tauw;
u1 = sqrt(u1);
plot(y2,u1,'r');hold on;
u1 = Fprof(1,:,WW).*Yprof(1,:,RHO)/tauw;
u1 = sqrt(u1);
plot(y2,u1,'g');hold on;



xlim([0 2])

MT = 0;
y3 = y2*delBL;
u1 = Yprof(1,:,U);
rho = Yprof(1,:,RHO);
u0 = u1(end);
uM = u1.*rho/u0/rhoin.*(1-u1/u0);
MT = trapz(y3,uM);

Re_theta = u0*rhoin*MT/mu_0

uM = 1-u1.*rho/u0/rhoin;
DT = trapz(y3,uM);
Re_D = u0*rhoin*DT/mu_0
%for i=1:ny/2
%    MT = MT + u1(i)/u0*(1-u1(i)/u0);
%end

 
figure(1);
xlabel('y+');
ylabel('u+');

%pause 

%[h,ypm] = hdrload('yplus_med.dat');    
%figure(1);hold on;
%plot(ypm(:,1),ypm(:,2),'r');


figure(2);
xlabel('y/delta');
ylabel('uu');

pprof(:,1) = Yprof(1,:,Y);
pprof(:,2) = y1;
pprof(:,3) = del + 0*y1;
pprof(:,4) = Yprof(1,:,U);
pprof(:,5) = u1;
pprof(:,6) = uvd;
pprof(:,7) = utau + 0*u1;
pprof(:,8) = Yprof(1,:,UU);
pprof(:,9) = Yprof(1,:,VV);
pprof(:,10) = Yprof(1,:,WW);
pprof(:,11) = delBL + 0*y1;
pprof(:,12) = tauw + 0*y1;

key1 = '%% < 1-5  > y(cm), y+, del, u(cm/s),u+';
key2 = '%% < 6-10 > uvd, utau, uu, vv, ww';
key3 = '%% < 11-12> delBL, tauw';
%dlmwrite(ofile,key1,'delimiter',' ');
%dlmwrite(ofile,key2 ,'delimiter',' ','-append');
%dlmwrite(ofile,pprof,'delimiter',' ','-append');

