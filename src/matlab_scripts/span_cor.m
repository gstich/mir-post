function span_cor


data = '/p/lscratchd/olson45/nozzle/';
gdir = [data,'post_proccoarse3d'];
gxy = load_planes_grid(gdir);

dir = [data,'post_proccoarse3d/'];
ifile = '~/miranda/chaos_4_x86_64_ib/bin/corr.cor';
var = 'u';
dz = .021;
nz = 128;
Lt = 1.78;
intv = [1200,1248];


ii = 0;
npts = 8*7*3*3;
corr(1:128,1:npts) = 0;
for i=intv(1):intv(2)
    ii = ii + 1;
    file = [ dir ,'post',num2str(i,'%04i'),'/',var,'-corr.dat' ]
    
    data1 = load(file);
    mm = 1;
    data = zeros(nz,npts);
    for p=1:npts
        for k=1:nz
            data(k,p) = data1(mm);
            mm = mm + 1;
        end
    end
    
    for j=1:size(data,2)
        corr(:,j) = corr(:,j) + get_corr(data(:,j));
    end
    
end

corr = corr / (intv(2)-intv(1)+1);

%figure(2);
%plot(corr(1:end/2,:))

Cor = corr(1:end/2,:);
%p_m = (sign(Cor) + 1)/2;
%Cor = Cor.*p_m;

for j=1:size(data,2)
    %taylor(j) = trapz(abs(corr(1:end/2,j)));
    ddf = -5/2*Cor(1,j) + 8/3*Cor(2,j) - 2/12*Cor(3,j);
    taylor(j)  = ( -.5*ddf )^(-.5);
    %taylor(j) = trapz(Cor(:,j));
end

Pmap = load(ifile);

np = Pmap(1,1);
mx = Pmap(1,2);
my = np/mx;
ii = 1;
for j=1:my
    for i=1:mx
        ii = ii+1;
        xmap(i,j) = Pmap(ii,1);
        ymap(i,j) = Pmap(ii,2);
        tay(i,j) = taylor(ii-1);
    end
end

for i=1:mx
    for j=1:my
        XX(i,j) = gxy( xmap(i,j), ymap(i,j),1);
        YY(i,j) = gxy( xmap(i,j), ymap(i,j),2);
    end
end

figure(2);
tay = tay*dz/Lt;
contourf(XX,YY,tay,32,'edgecolor','none');
colorbar;

%save '../data/span_corr/cor.mat' XX YY tay Pmap Cor

end



function corr = get_corr(data)


ave = mean(data);
data = data - ave;
nz = max(size(data));

% Periodic copy
dtmp(1:nz) = data(1:nz);
dtmp(nz+1:2*nz) = data(1:nz);

corr = zeros(nz,1);
for i=1:nz % Loop over lengths
    
    for j=1:nz
        
        corr(i) = corr(i) + dtmp(j)*dtmp(j+i-1);
        
    end
    
end

corr = corr / corr(1);

end