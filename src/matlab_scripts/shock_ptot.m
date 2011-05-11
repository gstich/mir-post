function shock_ptot


t1 = 350;
t2 = 665;
path = '/p/lscratchd/olson45/nozzle/post_procmedium3d';

ix1 = 400;
ixn = 600;
iy1 = 28;
iyn = 228;
ix = 500;
iy = 100;
var = 12;



for tt=t1:t2
    plane = load_planes_mir(path,'post',tt);

    Pdata(tt-t1+1) = plane(ix,iy,12);

    data = plane(ix1:ixn,iy1:iyn,12);
    tdata(:,:,tt-t1+1) = coarsen( data , 2);
end
   

plot(Pdata);
save pdata;


load('pdata');
load('Xs');
nt = size(xs,1);

%xs = xs - mean(xs);
Pdata = Pdata - mean(Pdata);

AVE = mean(tdata,3);

for i=1:size(tdata,3)
tdata(:,:,i) = tdata(:,:,i) - AVE;
end

%% Get the xs-p' correlation in tau > 0
count = zeros(nt,1);
tpcor = zeros(size(tdata));
ptcor = zeros(size(tdata));
for t=1:nt

    for i=1:nt
        
        off = t+i-1;
        if (off <= nt) 
            tpcor(:,:,t) = xs(t) * tdata(:,:,off);
            ptcor(:,:,t) = tdata(:,:,t) * xs(off);
            count(t) = count(t) + 1;
        end
        
    end
    tpcor(:,:,t) = tpcor(:,:,t) / count(t);
    ptcor(:,:,t) = ptcor(:,:,t) / count(t);
    
end

corr(:,:,1:nt) = ptcor(:,:,nt:-1:1);
corr(:,:,nt+1:2*nt) = tpcor(:,:,1:nt);
dd(1:2*nt) = corr(80,50,:);
plot(dd)

end