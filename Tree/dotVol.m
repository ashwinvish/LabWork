function [area] = dotVol( vol1, vol2 , soma1, soma2, res )
%dotVol dotprodict of the interaction of two volumes
%   vol1, soma1 presynaptic cell volume, soma
%   vol2, soma2 postsynaptic cell volume, soma
%   res is the downsampling factor

% downsampling by res
soma1 = soma1./res;
soma2 = soma2./res; 

tempxy = dot(vol1,vol2,3);
tempxz = dot(vol1,vol2,2);
tempyz = dot(vol1,vol2,1);

%cmap = cbrewer('div','Spectral', size(vol1,1));

if ~(isempty(find((tempxy>0))& isempty(find(tempxz>0))& isempty(find(tempyz>0)))>0)
    area = threeView([],soma1,jet,tempxy,tempxz,tempyz, soma2);
    area = area*res*res; % back to scale, since volumes are calculated at a downsampled resolution
else
    area = 0;
end

end

