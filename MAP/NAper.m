function [na] = NAper(Acc , PixSize)
%   NAper, Numerical Aperture of SEM given Accelerating 
%   voltage and Pixel size

na = 0.752 / (PixSize* (Acc*1000)^0.5);


end

