function wall_spectra


res = 'coarse';
%res = 'medium';

switch(res)
    case('coarse');
        load('DATA/Pw_cor.mat');
        dx = 1;
        x1 = 100;
        xn = 460;
        xs = 5;
        %Ptop = P;
    case('medium')
        load('DATA/Pw_med.mat');
        dx = 1;
        x1 = 170;
        xn = 660;
        xs = 8;
        %Ptop = Pbot;
end



zs = 1;

L = size(Pw,3);
NFFT = 2^nextpow2(L);
Nt = NFFT/2+1;
pgrid = zeros( (xn-x1)/xs, Nt);
ptime = zeros( L,1);
pspec = zeros(Nt,1);
xcount = 0;
for i=x1:xs:xn
    

    pspec = 0*pspec;
    zcount = 0;
    for k=zs:zs:size(Pw,2)
        
        ptime(1:end) = Pw(i,k,:);
        
        a = get_spectra(ptime);
        %a = a.*a;
        
        pspec = pspec + a(2,:)';
        zcount = zcount + 1;
    end

    xcount = xcount + 1;
    pgrid(xcount,1:end) = pspec(1:end) / zcount;
    
end

f = a(1,:);

% Compensated spectra
Ht = 1.78;
Up = 32940*1.603;

FF = zeros(xcount,Nt);
XX = FF;
for i=1:xcount
    pgrid(i,:) = pgrid(i,:).*f;
    Pnorm = sum(pgrid(i,:));
    pgrid(i,:) = pgrid(i,:) / Pnorm;
    FF(i,:) = f; %log(f);
    XX(i,:) = (i-1)*dx; 
end
pgrid = pgrid/Nt;
% non- by max val
pgrid = pgrid / max(max(pgrid));

xcount = 0;
for i=x1:xs:xn
    ii = xcount + 1;
    XX(ii,:) = X(i)/Ht;
    xcount = xcount + 1;
end
%hold all;
%for i=1:xcount
%    plot(FF(i,:),pgrid(i,:));
%end

%for i=1:xcount
%    % Filter this term
%    b = pgrid;
%    for k=1:1
%        for ii=2:max(size(b))-1
%            b(:,ii)=(b(:,ii-1)+b(:,ii+1))/2;
%        end
%    end
%    pgrid = b;
%end

write_visit(['DATA/sim_',res,'_spec.tec'],XX(:,2:end),FF(:,2:end),pgrid(:,2:end));
figure(2);
contourf(XX,FF,pgrid,16,'edgecolor','none');
xray = flipud(gray);
colormap(xray);

%save WS.mat




end


function a = get_spectra(ptot)

Ht = 1.78;
Up = 32940*1.603;
L=size(ptot,1);
Fs = 200e3;
hw = hanning(L,'periodic');

y = ptot - mean(ptot);
y = y/10;   % Convert to Pascals
y = y;%.*y;
y = y.*hw;


NFFT = 2^nextpow2(L); % Next power of 2 from length of y
%figure(2)
Y = (fft(y,NFFT))/L;
f = Fs/2*linspace(0,1,NFFT/2+1);

%% Divide by frequency

a(1,:) = f*Ht/Up;
a(2,:) = 2*abs(Y(1:NFFT/2+1));


end