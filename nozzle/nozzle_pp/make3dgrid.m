clear all;
clc;

% Make a 3d grid using 2 planar grids
nx = 256;
ny = 64;
nz = 64;

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
for j=1:nz
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

%Nb = 500;
%yblnd = linspace(0,yH,Nb);
%for i=1:Nb
%    if (yblnd(i) <= yL)
%        blend(i,1) = 0;
%    else
%        % Try linear first
%        
%       d = yblnd(i)-yL;
%        L = yH-yL;
%        A = 1/L^2;
%        blend(i,1) = d/L;     % Linear
%        blend(i,1) = A*d^2;   % Parabolic
%    end
%end

%Nfil = 10;
%for i=1:Nfil
%   blend = gfilter(blend);
%end
figure(12);hold all;
xline = gXZx(:,1); % Make xline dependent on y-location
for i=1:nx
    for j=1:ny
        
        yht = abs(gXYy(i,j));
        %if (yht<=yL)
        %    xline = gXZx(:,1); % Make xline dependent on y-location
        %else
            % Try linear first
            d = yht-yL;
            L = yH-yL;
            A = 1/L^2;
            L1 = .707*yH-yL;
            mid = .707*yH-L1/2;
            %blend = d/L;   % Linear
            %blend = A*d^2; % Parabolic
            blend = .5*(1+tanh((yht-mid)/(L1/4)));
            xline = gXZx(:,1)*(1-blend) + gXZ2x(:,1)*blend;
        %end
        
        %if(i>nx-5)
        %figure(12);
        %plot(yht,blend,'bo');
        %end
        % Make very smooth transition (tanh)
        %smth = .5*(1+tanh((yht-yL)/yL*(1/5)));
        
        [a,b] = min(abs(gXYx(i,j) - xline));
        if (b==1)
            b = 2;
        elseif (b==nx)
            b = nx-1;
        end
        vm = ( gXYx(i,j)-xline(b-1) );
        v0 = ( gXYx(i,j)-xline(b) );
        vp = ( gXYx(i,j)-xline(b+1) );
        
        if ( sign(vm) == -sign(v0))
            l1=abs(vm);
            l2=abs(v0);
            L = l1+l2;
            pt1=b-1;
            pt2=b;
            xit = (l2*pt1+l1*pt2)/L;
        elseif( sign(v0) == -sign(vp) )
            l1=abs(v0);
            l2=abs(vp);
            L = l1+l2;
            pt1=b;
            pt2=b+1;
            xit = (l2*pt1+l1*pt2)/L;
        elseif (sign(vm)==0)
            xit = b-1;
        elseif (sign(v0)==0)
            xit = b;
        elseif (sign(v0)==0)
            xit = b+1;
        end
   
        xmap(i,j) = xit;
        
    end
end
      
%pause;

gXZNx = gXZ2x*0;
gXZNz = gXZ2x*0;

for i=1:nx
    i
    for j=1:ny
    for k=1:nz
        
        % Project in -y direction and see what mapped x-point where at
        xx = xmap(i,j);
        x1 = floor(xx);
        x2 = ceil(xx);
        w2 = rem(xx,x1);
        
        yht = abs(gXYy(i,j));
        if (yht<=yL)
            gXZNx = gXZx; 
            gXZNz = gXZz; 
        else
            % Try linear first
            d = yht-yL;
            L = yH-yL;
            gXZNx = gXZx*(1-d/L) + gXZ2x*d/L;
            gXZNz = gXZz*(1-d/L) + gXZ2z*d/L;
        end
        
        
        if (x1==x2)
            xoff = gXZNx(x1,k) - gXZx(x1,1);
            zoff = gXZNz(x1,k) - gXZz(x1,1);
        else
            x1off = gXZNx(x1,k) - gXZNx(x1,1);
            z1off = gXZNz(x1,k) - gXZNz(x1,1);
        
            x2off = gXZNx(x2,k) - gXZNx(x2,1);
            z2off = gXZNz(x2,k) - gXZNz(x2,1);
            
            xoff = (1-w2)*x1off + w2*x2off;
            zoff = (1-w2)*z1off + w2*z2off;
        end
        
        
        X(i,j,k) = gXYx(i,j) + xoff;
        Z(i,j,k) = zoff;
        Y(i,j,k) = gXYy(i,j);
        
    end
    end
end
  

figure(4 );
for kk=1:1
mesh(X(:,:,nz+1-kk),Y(:,:,nz+1-kk),Z(:,:,nz+1-kk),'edgecolor','black');hold on;
end

view(2)
axis equal


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





