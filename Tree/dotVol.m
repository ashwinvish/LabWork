function [area ] = dotVol( vol1, vol2 , soma1, soma2, res )
%dotVol dotprodict of the interaction of two volumes
%   vol1, soma1 presynaptic cell volume, soma
%   vol2, soma2 postsynaptic cell volume, soma
%   res is the downsampling factor

Soma = [soma1./res;soma2./res]; % downsampling by res

tempxy = dot(vol1,vol2,3);
tempxz = dot(vol1,vol2,2);
tempyz = dot(vol1,vol2,1);


if ~(isempty(find((tempxy>0))& isempty(find(tempxz>0))& isempty(find(tempyz>0)))>0)
    area = threeView([],Soma,jet,tempxy,tempxz,tempyz);
    area = area*10e6;
else
    area = 0;
end

end

