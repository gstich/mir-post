function [xs,ps] = get_shock_location(x,press)

    val = press;
    for j=1:50
        val = gfilter(val);
    end
    
    
    %% Try fitting parabola to min.
    off = 1;
    [Pmin,Nmin] = min(val);
    X1 = x(Nmin);
    Y1 = Pmin;
    X2 = x(Nmin - off);
    Y2 = val(Nmin - off);
    X3 = x(Nmin + off);
    Y3 = val(Nmin + off);
    
    DD = X1*(X3^2-X2^2)-X2*X3^2 + X2^2*X3 + X1^2*(X2-X3);
    A = X1*(Y3-Y2)-X2*Y3+X3*Y2+(X2-X3)*Y1;
    A = A/DD;
    B = X1^2*(Y3-Y2)-X2^2*Y3+X3^2*Y2+(X2^2-X3^2)*Y1;
    B = B/DD;
    
    xs =  B / (2*A);
    ps = Pmin;

end