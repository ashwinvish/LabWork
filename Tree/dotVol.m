function [area ] = dotVol( vol1, vol2 ,soma1,soma2,res)
%dotVol dotprodict of the interaction of two volumes
%   vol1, soma1 presynaptic cell volume, soma
%   vol2, soma2 postsynaptic cell volume, soma
Soma = [soma1./res;soma2./res];

tempxy = dot(vol1,vol2,3);
tempxz = dot(vol1,vol2,2);
tempyz = dot(vol1,vol2,1);

if ~(isempty(find((tempxy>0))& isempty(find(tempxz>0))& isempty(find(tempyz>0)))>0)
    area = threeView([],jet,tempxy,tempxz,tempyz);
    %cell1 - red
    %plot XZ soma location
    plot(160-(Soma(1,1)), 80 -(Soma(1,3)),'Marker','o', 'MarkerFaceColor',rgb('gray'), 'MarkerEdgeColor','r');
    %plot XY soma location
    plot(160-(Soma(1,1)), 90 + (Soma(1,2)),'Marker','o', 'MarkerFaceColor',rgb('gray'), 'MarkerEdgeColor','r');
    %plot YZ soma location
    plot(160+(Soma(1,3)), 90 + (Soma(1,2)),'Marker','o', 'MarkerFaceColor',rgb('gray'), 'MarkerEdgeColor','r');
    % cell2 - green
    %plot XZ soma location
    plot(160-(Soma(2,1)), 80 -(Soma(2,3)),'Marker','o', 'MarkerFaceColor','w', 'MarkerEdgeColor','g');
    %plot XY soma location
    plot(160-(Soma(2,1)), 90 + (Soma(2,2)),'Marker','o', 'MarkerFaceColor','w', 'MarkerEdgeColor','g');
    %plot YZ soma location
    plot(160+(Soma(2,3)), 90 + (Soma(2,2)),'Marker','o', 'MarkerFaceColor','w', 'MarkerEdgeColor','g');
    axis off;
    axis vis3d;
else
    area = 0;
end

end

