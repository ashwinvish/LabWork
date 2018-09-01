function [newMap] = NLmap(map,dataMax, dataMin, centerPoint, scalingIntensity)

cMap = map; % 256 colors;

%Then perform some operations to create your colormap. 
%I have done this by altering the indices “x” at which each existing color lives, 
%and then interpolating to expand or shrink certain areas of the spectrum.

x = 1:length(cMap); 
x = x - (centerPoint-dataMin)*length(x)/(dataMax-dataMin);
x = scalingIntensity * x/max(abs(x));

%Next, select some function or operations to transform the original linear indices into nonlinear. 
%In the last line, I then use “interp1” to create the new colormap from the original colormap and
%the transformed indices.

x = sign(x).* exp(abs(x));
%x = log(abs(x));
x = x - min(x); 
x = x*511/max(x)+1; 
newMap = interp1(x, cMap, 1:512);
end

