function dF = ddx(F)

n = max (size(F) );


    I = [2:n-1];
    Ip = I + 1;
    Im = I -1;

    dF(I,:) = 1/2*( F(Ip,:) - F(Im,:) );
    
    dF(1,:) =-F(1,:) + F(2,:) ; 
    dF(n,:) =    F(n,:) -F(n-1,:);
    

    



end