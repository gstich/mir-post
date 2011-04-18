
clc;

% Some parameters for ease of use
X=1; Y=2; U=3; V=4; W=5; RHO=6; P=7; T=8;


t1 = 600;
t2 = 602;
Nwall = 80;
buff = 5;
path = '/p/lscratchd/olson45/nozzle/nozzlecoarse3d';

if  ~exist('Pmean')
    Pmean = sum_planes(path,'mean',t1,t2);
end
if  ~exist('fluc')
    fluc = sum_planes(path,'flc' ,t1,t2);
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
Mach = 1.0;
ReBL = 5e3;
delBL= 1.5e-1;
gamma = 1.4;
Pin  =   Pmean(1,ny/2, P );
rhoin  = Pmean(1,ny/2,RHO);
Uin = Mach*sqrt(Pin*gamma/rhoin);
mu_0 = Uin * rhoin * delBL / ReBL;




% Sutherland's law
T_0 = 273.15;  % Kelvin reference temperature
ST = 110.4;    % Sutherland temperature

T_wall = Yprof(1,1,T);
mu_w = mu_0 * ((T_wall/T_0)^(3.0D0/2.0D0) * (T_0 + ST) / (T_wall + ST) );
mu_w = mu_w;% * 1.8;

% Scaling
rho_w = Yprof(1,1,RHO);
dudy = ( Yprof(1,2,U) - Yprof(1,1,U) )/ ( Yprof(1,2,Y) - Yprof(1,1,Y) );
tauw = mu_w * dudy;
utau = sqrt( tauw / rho_w );
del = mu_w / ( rho_w * utau );

figure(1);%clf;
y1 = Yprof(1,:,Y)-Yprof(1,1,Y);
u1 = Yprof(1,:,U);
y1 = y1/del;
u1 = u1/utau;
semilogx( y1,u1,'bo');
x1 = logspace(-1,1.2,20);
x2 = logspace(.8,3,20);
k = .42;
C = 5.2;
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



figure(2);
y1 = y1*del/delBL;
u1 = Fprof(1,:,U).*Yprof(1,:,RHO)/tauw;
u1 = sqrt(u1);
plot(y1,u1,'b');hold on;
u1 = Fprof(1,:,V).*Yprof(1,:,RHO)/tauw;
u1 = sqrt(u1);
plot(y1,u1,'r');hold on;
u1 = Fprof(1,:,W).*Yprof(1,:,RHO)/tauw;
u1 = sqrt(u1);
plot(y1,u1,'g');hold on;



xlim([0 3])


