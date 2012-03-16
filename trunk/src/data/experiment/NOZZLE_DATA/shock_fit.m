function xs = shock_fit(p,x,ps1,ps2,xg)

% Compute the location of the shock for a given pressure profile by the
% method as outlined in Johnson and Papamoschou, 2010.

del = (ps2 - ps1) / max( (p(2:end) - p(1:end-1))/(x(2:end) - x(1:end-1)));

lambda = 0.05;

xs = xg;
res = 1;
Eold = 0;
while res > 1e-6

    [E,e] = error(xs,x,p,ps1,ps2,del);
    
    J = dPfun(x,xs,del,ps1,ps2);
    dxs = J*e';
    xs = xs + lambda*dxs;
    if(xs < x(1))
        p(1) = ps1;
        if (xs < .8*x(1))
            xs = .8*x(1);
        end
    end
    
    if(xs > x(end))
        p(end) = ps2;
        if (xs>1.2*x(end))
            xs = 1.2*x(end);
        end
    end
    
    
    res = abs( E - Eold );
    Eold = E;
end
    
end



function [E,e] = error(xs,x,p,ps1,ps2,del)

prof = Pfun(x,xs,del,ps1,ps2);
for i=1:4
    e(i) = p(i) - prof(i);
end

E = sum( e.*e );

end


function p = Pfun(x,xs,del,ps1,ps2)

p = (ps2+ps1)/2 + (ps2-ps1)/2 * tanh( 2/del*(x-xs) );

end


function J = dPfun(x,xs,del,ps1,ps2)

J = - (ps2-ps1)/del * sech( 2*(x-xs)/del).^2;

end