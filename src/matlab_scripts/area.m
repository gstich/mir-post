clc;
clear all;
close all;

% Some parameters for ease of use
X=20; Y=21; U=1; V=2; W=3; RHO=5; P=4; T=6;UU=7;VV=8;WW=9;
MU=16;MUa=17;MUb=18;MUc=19;


%% Parameters
res = 'coarse';
%res = 'medium';

P0 = 1.0e6;
Ht = 1.78;
switch (res)
    case 'coarse'
        path = '/p/lscratchd/olson45/nozzle/post_proccoarse3d';
        nx = 512;
        ny = 128;
        t1 = 1000;  
        tf = 2548;
        xbound = [0,11];
        
    case 'medium'
        path = '/p/lscratchd/olson45/nozzle/post_procmedium3d';
        nx = 768;
        ny = 256;
        t1 = 400;  
        tf = 1100;
        xbound = [0,11];
end


data = sum_planes(path,'post',t1,t1);
tmp = data(2:end,:,X);
[a,b] = min(abs(tmp(:,1)-xbound(1)));
ixbound(1) = b;
[a,b] = min(abs(tmp(:,1)-xbound(2)));
ixbound(2) = b;
x_c = tmp(ixbound(1):ixbound(2),:);
y_c = data(ixbound(1):ixbound(2),:,Y);

wall = y_c(:,end)*2;
clear tmp;
tmp = data;clear data;
data = tmp(ixbound(1):ixbound(2),:,:);

 


Nxx = size(data,1);
map = zeros(Nxx,size(data,1));
weight = map;
clear data;
for tt=t1:tf
  it = tt-t1+1;
  tmp = load_planes_mir(path,'post',tt);
  data = tmp(ixbound(1):ixbound(2),:,:);

  for i=1:size(data,1)
  
    u1 = data(i,:,U)';
    u2 = 1/2*(1+sign(u1));
    u1 = u2 / max(max(u2));
    u2 = u1;%.^4;
    
    W(i) = trapz(y_c(i,:)',u2);
  end
 
  
  
  
  % Get the min pressure
  p = data(:,ny/2,P);
  [pmin,imin] = min(p);
  xs = x_c(imin,ny/2);
  [xs,pmin] = get_shock_location(x_c(:,ny/2),p);
  
  % Add to the map... bins
  off = (wall'-W)/Ht;
  %map(imin,:) = map(imin,:) + off;
  %weight(imin,:) = weight(imin,:) + 1;
  
  width = Ht/10;   % width of the gauss
  gauss = zeros(size(data,1),size(data,1));
  wght = gauss;
  for i=1:size(data,1)
      wt = exp( -( xs - x_c(i,ny/2))^2 / width^2 );
      gauss(i,:) = wt * off;
      wght(i,:) = wt;
  end
  
  % Add to map ... gauss'
  map = map + gauss;
  weight = weight + wght;
  
  
  %figure(1);clf;
  %off = (wall'-W)/Ht;
  %plot(x_c(:,ny/2),off);hold on;
  %plot(xs,off(imin),'ro');
  %ylim([0 .5])
end

map = map./weight;

%figure(2);
%contour(x_c(:,ny/2),x_c(:,ny/2),map');
%ylim([6 11.5]);xlabel('Shock Location')
%xlim([6.5 10]);ylabel('Nozzle Length')

figure(3);
contourf(x_c(:,ny/2),x_c(:,ny/2),map',32,'edgecolor','none');
ylim([6 11.5]);xlabel('Shock Location')
xlim([6.5 10]);ylabel('Nozzle Length')

figure(3);hold on;
plot(x_c(:,ny/2),weight(:,1)/max(max(weight))+6.5,'k-','Linewidth',2)



