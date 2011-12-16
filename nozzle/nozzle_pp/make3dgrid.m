clear all;
clc;

% Make a 3d grid using 2 planar grids
nx = 256;
ny = 64;
nz = 64;

%itype = 'nearest';
itype = 'spline';
%itype = 'linear';

gXY = load('nozzleXY.grid');

cc = 1;
for j=1:ny
    for i=1:nx
        gXYx(i,j) = gXY(cc,1);
        gXYy(i,j) = gXY(cc,2);
        cc = cc + 1;
    end
end

figure(1);
mesh(gXYx,gXYy,ones(256,64),'edgecolor','black')
view(2)
axis equal



gXZ = load('nozzleXZ.grid');

cc = 1;
for j=nz:-1:1
    for i=1:nx
        gXZx(i,j) = gXZ(cc,1);
        gXZz(i,j) = gXZ(cc,2);
        cc = cc + 1;
    end
end

figure(2);
mesh(gXZx,gXZz,ones(256,64),'edgecolor','black')
view(2)
axis equal


gXZ2 = load('nozzleXZ2.grid');

cc = 1;
for j=1:nz
    for i=1:nx
        gXZ2x(i,j) = gXZ2(cc,1);
        gXZ2z(i,j) = gXZ2(cc,2);
        cc = cc + 1;
    end
end

figure(3);
mesh(gXZ2x,gXZ2z,ones(256,64),'edgecolor','black')
view(2)
axis equal



% Loop over planes of constant K (zdir) get the nominal offset at y==0 and
% apply to all
X = zeros(nx,ny,nz);
Y = zeros(nx,ny,nz);
Z = zeros(nx,ny,nz);

yL = gXYy(1,end) * 4.8;
yH = gXYy(end,end);


figure(12);hold all;
xline = gXZx(:,1); % Make xline dependent on y-location
for i=1:nx
    for j=1:ny
        
        yht = abs(gXYy(i,j));
        blend = blendY(yht,yL,yH);
        xline = gXZx(:,1)*(1-blend) + gXZ2x(:,1)*blend;
        
        [a,b] = min(abs(gXYx(i,j) - xline));
        if (b==1)
            b = 2;
        elseif (b==nx)
            b = nx-1;
        end
        vm = ( gXYx(i,j)-xline(b-1) );
        v0 = ( gXYx(i,j)-xline(b) );
        vp = ( gXYx(i,j)-xline(b+1) ); 
        
        cc = Finterp([vm,v0,vp],[b-1,b,b+1],0.0,itype);

        xmap(i,j).cc = cc;
        xmap(i,j).b = b;
        
    end
end
      

gXZNx = gXZ2x*0;
gXZNz = gXZ2x*0;

for i=1:nx
    i
    for j=1:ny
    for k=1:nz
        
        xx = xmap(i,j).b;
        cc = xmap(i,j).cc;
        xm = xx - 1;
        xp = xx + 1;
        
        
        yht = abs(gXYy(i,j));
        blend = blendY(yht,yL,yH);

        
        xrange = [xm,xx,xp];
        gXZNx(xrange,k) = gXZx(xrange,k)*(1-blend) + gXZ2x(xrange,k)*blend;
        gXZNz(xrange,k) = gXZz(xrange,k)*(1-blend) + gXZ2z(xrange,k)*blend;    
        
        xmoff = gXZNx(xm,k) - gXZNx(xm,1);
        zmoff = gXZNz(xm,k) - gXZNz(xm,1);
        
        xxoff = gXZNx(xx,k) - gXZNx(xx,1);
        zzoff = gXZNz(xx,k) - gXZNz(xx,1);
        
        xpoff = gXZNx(xp,k) - gXZNx(xp,1);
        zpoff = gXZNz(xp,k) - gXZNz(xp,1);
        
        xoff = Finterp( [xm,xx,xp],[xmoff,xxoff,xpoff],cc,itype);
        zoff = Finterp( [xm,xx,xp],[zmoff,zzoff,zpoff],cc,itype);

        
        X(i,j,k) = gXYx(i,j) + xoff;
        Z(i,j,k) = zoff;
        Y(i,j,k) = gXYy(i,j);
        
    end
    end
end
  

if (1==1)
        
% Write tecplot file
tmp = zeros(256*64*64,4);
cc = 1;
for k=1:nz
    k
    for j=1:ny
        %for i=1:nx
            tmp(cc:cc+nx-1,1:4) = [X(1:nx,j,k),Y(1:nx,j,k),Z(1:nx,j,k),ones(nx,1)];
            cc = cc + nx;
        %end
    end
end

vizname = 'visit2.tec';

dlmwrite(vizname,'VARIABLES="X","Y","Z","Pressure" ZONE I=','delimiter','')
dlmwrite(vizname,size(X,1),'delimiter','','-append')
dlmwrite(vizname,',J=','delimiter','','-append')
dlmwrite(vizname,size(X,2),'delimiter','','-append')
dlmwrite(vizname,',K=','delimiter','','-append')
dlmwrite(vizname,size(X,3),'delimiter','','-append')
dlmwrite(vizname,',F=POINT','delimiter','','-append')

dlmwrite(vizname,tmp,'delimiter',' ','-append');

end





