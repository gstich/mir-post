function weight = blendY(y,y1,y2)

               
d = y-y1;
L = y2-y1;          
A = 1/L^2;          
L1 = .707*y2-y1;          
mid = .707*y2-L1/2;  

if (y<=y1)
    blend = 0;
else
    blend = d/L;   % Linear
    %blend = A*d^2; % Parabolic            
end

blend = .5*(1+tanh((y-mid)/(L1/4)));  % Tanh

weight = blend;


end