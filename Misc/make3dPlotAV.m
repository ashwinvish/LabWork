function make3dPlotAV(s, cmap)
% make3dPlotAV Summary of this function goes here
%   Detailed explanation goes here
[y x z] = size(s)/2;

Ix = toColormappedChannelAV(squeeze(s(y,:,:)), cmap, 0.1); imshow(Ix);
Iy = toColormappedChannelAV(squeeze(s(:,x,:)), cmap, 0.1); imshow(Iy);
Iz = toColormappedChannelAV(squeeze(s(:,:,z)), cmap, 0.1); imshow(Iz);

do3Dplot(Ix,Iy,Iz);
view([90 80 40]);
light('Position',[5 8 12],'Style','local')
axis off
end

