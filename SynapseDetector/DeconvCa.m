function [dstaf] = DeconvCa(staf)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

delt = 0.05;
tau_cirf = 1.9;

dlen = floor(5/delt);           % length of kernel (in index units)
Tnew = linspace(0,dlen*delt,dlen);
ker = exp(-(Tnew-Tnew(1))/tau_cirf);
zpad = zeros(dlen-1,1);         % zero-padding


dstaf = (tau_cirf/delt)*deconv([staf;zpad],ker);
dstaf = medfilt1(dstaf,20);

end

