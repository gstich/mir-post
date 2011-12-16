function yi = Finterp(xx,yy,xi,order)


switch order
    case 'nearest';
        
        [a,b] = min( abs(xx-xi) );
        yi = yy(b);
        
    case 'linear';
        
        if ( xi < xx(2) )
            yi = yy(2) + ( yy(3)-yy(2) )/( xx(3)-xx(2) )*(xi-xx(2));
        else
            yi = yy(1) + ( yy(2)-yy(1) )/( xx(2)-xx(1) )*(xi-xx(1));
        end
            
    case 'spline';
        
        a = (xx(1)*(yy(3)-yy(2))-xx(2)*yy(3)+xx(3)*yy(2)+(xx(2)-xx(3)) ... 
            *yy(1))/(xx(1)*(xx(3)^2-xx(2)^2)-xx(2)*xx(3) ...
            ^2+xx(2)^2*xx(3)+xx(1)^2*(xx(2)-xx(3)));
        b = -(xx(1)^2*(yy(3)-yy(2))-xx(2)^2*yy(3)+xx(3)^2*yy(2) ...
            +(xx(2)^2-xx(3)^2)*yy(1))/(xx(1)*(xx(3)^2-xx(2)^2)-xx(2) ...
            *xx(3)^2+xx(2)^2*xx(3)+xx(1)^2*(xx(2)-xx(3)));
        c = (xx(1)*(xx(3)^2*yy(2)-xx(2)^2*yy(3))+xx(1)^2*(xx(2)*yy(3)- ...
            xx(3)*yy(2))+(xx(2)^2*xx(3)-xx(2)*xx(3)^2)*yy(1))/(xx(1)* ...
            (xx(3)^2-xx(2)^2)-xx(2)*xx(3)^2+xx(2)^2*xx(3)+ ... 
            xx(1)^2*(xx(2)-xx(3)));
        
        yi = a*xi^2 + b*xi + c;
        
end



end