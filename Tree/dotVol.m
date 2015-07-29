function [ tempxy,tempxz,tempyz ] = dotVol( vol1, vol2 ,soma1,soma2,res)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
Soma = [soma1./res;soma2./res];

tempxy = dot(vol1,vol2,3);
tempxz = dot(vol1,vol2,2);
tempyz = dot(vol1,vol2,1);
threeView([],jet,tempxy,tempxz,tempyz); 

%cell1
%plot XZ soma location
plot(160-(Soma(1,1)), 80 -(Soma(1,3)),'Marker','o', 'MarkerFaceColor','w', 'MarkerEdgeColor','r');
%plot XY soma location
plot(160-(Soma(1,1)), 90 + (Soma(1,2)),'Marker','o', 'MarkerFaceColor','w', 'MarkerEdgeColor','r');
%plot YZ soma location
plot(160+(Soma(1,3)), 90 + (Soma(1,2)),'Marker','o', 'MarkerFaceColor','w', 'MarkerEdgeColor','r');

% cell2
%plot XZ soma location
plot(160-(Soma(2,1)), 80 -(Soma(2,3)),'Marker','o', 'MarkerFaceColor','w', 'MarkerEdgeColor','g');
%plot XY soma location
plot(160-(Soma(2,1)), 90 + (Soma(2,2)),'Marker','o', 'MarkerFaceColor','w', 'MarkerEdgeColor','g');
%plot YZ soma location
plot(160+(Soma(2,3)), 90 + (Soma(2,2)),'Marker','o', 'MarkerFaceColor','w', 'MarkerEdgeColor','g');

axis off;
axis vis3d;
end

