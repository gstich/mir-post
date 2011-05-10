function shock_ptot


t1 = 500;
t2 = 620;
path = '/p/lscratchd/olson45/nozzle/post_procmedium3d';

ix = 500;
iy = 100;
var = 12;

for tt=t1:t2
    plane = load_planes_mir(path,'post',tt);

    Pdata(tt-t1+1) = plane(ix,iy,12);

end
   

plot(Pdata);
save pdata;

load('Xs');
nt = size(xs,2);

xs = xs - mean(xs);
Pdata = Pdata - mean(Pdata);

%% Get the xs-p' correlation in tau > 0
count = zeros(nt,1);
tpcor = zeros(nt,1);
ptcor = zeros(nt,1);
for t=1:nt

    for i=1:nt
        
        off = t+i-1;
        if (off <= nt) 
            tpcor(t) = xs(t) * Pdata(off);
            ptcor(t) = Pdata(t) * xs(off);
            count(t) = count(t) + 1;
        end
        
    end
end

corr = [ ptcor(end:-1:1) ; tpcor ];

end