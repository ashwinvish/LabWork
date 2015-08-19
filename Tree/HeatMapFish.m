function [plotVol,m,I] = HeatMapFish( ksize,res, PointCloud, Soma, CellID ,figure )
%HeatMap used to plot the heatmap of a PointCloud for a CellID 

%   PointCloud [X Y Z]  is the points of cloud over which the density is computed,
%   e.g clould be presynaptic coordinates.
%   ksize is the size of the kernel in same unit as coordinates ( here in
%   nm)
%   res is the downsampling factor
%   Soma (x y z) is the coordinates of the cell soma
%   CellID is the string to ID the cell
%   figure, true or false if figure is needed

% downsample the kernel size and the Soma coordinates
ksize = ksize/res;
Soma = Soma./res;

% set up 3D gaussian kernel
Gwin = gausswin(ksize+1);
K = Gwin*Gwin';

for i = 1:(ksize+1)
    K3D(:,:,i) = K(:,:)*Gwin(i);
end

% size of volume is prefixed
vol = zeros( ksize + (14e4-0)/res + ksize , ksize + (2.5e5-0)/res + ksize , ksize + (6e4-0)/res + ksize);

% pad with ksize/2 pixels to all the axes to avoid edge effects
for n = 1:size(PointCloud,1);
    vol(round(PointCloud(n,1)/res + ksize/2) + (-ksize/2:ksize/2), round(PointCloud(n,2)/res + ksize/2) + (-ksize/2:ksize/2), round(PointCloud(n,3)/res + ksize/2 +1)+ (-ksize/2:ksize/2)) = ...
        vol(round(PointCloud(n,1)/res + ksize/2) + (-ksize/2:ksize/2), round(PointCloud(n,2)/res + ksize/2) + (-ksize/2:ksize/2),round(PointCloud(n,3)/res + ksize/2 +1)+ (-ksize/2:ksize/2))  + K3D;
end

% normalize volume
vol = vol./max(vol(:));
% remove padding from the edges
plotVol = vol(20+ksize/2 : 140+ksize/2, 60+ksize/2 : 250+ksize/2, ksize/2:60+ksize/2);


if figure == true
    threeView(plotVol,Soma, jet);
    title(CellID);
    [m,I] = max(plotVol(:));
    
else
    [m,I] = max(plotVol(:));
end

end

