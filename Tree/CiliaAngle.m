function [ ag ] = CiliaAngle( x,y )
%CILIAANGLE plots the angle between the primary Cilia and the CellSomata
    z = (x + y.*sqrt(-1)).';
    [th,r] =  cart2pol(real(z),imag(z));
    ag = rad2deg(th);

end

