function [ vol,m,I] = HeatMapFish( ksize,res, PointCloud, Soma, CellID  )
%HeatMap used to plot the heatmap of a PointCloud for a CellID centered at
%Soma
%   PointCloud [X Y Z]  is the point of clouds over which the density is computed,
%   e.g clould be presynaptic coordinates.
%   Soma (x y z) is the coordinates of the cell soma
%   CellID is the string to ID the cell

Gwin = gausswin(ksize+1);
K = Gwin*Gwin';

for i = 1:(ksize+1)
    K3D(:,:,i) = K(:,:)*Gwin(i);
end

vol = zeros((15e4-0)/res,(2.6e5-0)/res,(8e4-0)/res); % size of volume is prefixed

% add 10000/res pixels to the zaxis to avoid edge effects
for n = 1:size(PointCloud,1);
    vol(round(PointCloud(n,1)/res) + (-ksize/2:ksize/2), round(PointCloud(n,2)/res) + (-ksize/2:ksize/2), round(PointCloud(n,3)/res + 10000/res)+ (-ksize/2:ksize/2)) = ...
        vol(round(PointCloud(n,1)/res) + (-ksize/2:ksize/2), round(PointCloud(n,2)/res) + (-ksize/2:ksize/2),round(PointCloud(n,3)/res + 10000/res)+ (-ksize/2:ksize/2))  + K3D;
end

threeView(vol,jet);
title(CellID);

%plot XZ soma location
plot(160-(Soma(1,1)/res), 80 -(Soma(1,3))/res,'Marker','o', 'MarkerFaceColor','w', 'MarkerEdgeColor','w');
%plot XY soma location
plot(160-(Soma(1,1)/res), 90 + (Soma(1,2))/res,'Marker','o', 'MarkerFaceColor','w', 'MarkerEdgeColor','w');
%plot YZ soma location
plot(160+(Soma(1,3)/res), 90 + (Soma(1,2))/res,'Marker','o', 'MarkerFaceColor','w', 'MarkerEdgeColor','w');

axis off
axis vis3d

[m,I] = max(vol(:));

end

