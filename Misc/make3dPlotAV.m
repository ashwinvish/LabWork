function make3dPlotAV(s)
% make3dPlotAV Summary of this function goes here
%   Detailed explanation goes here
coord = round(size(s)/2);
Ix = toColormappedChannelAV(squeeze(s(coord(2),:,:)), cmap, 0.1); imshow(Ix);
Iy = toColormappedChannelAV(squeeze(s(:,coord(1),:)), cmap, 0.1); imshow(Iy);
Iz = toColormappedChannelAV(squeeze(s(:,:,coord(3))), cmap, 0.1); imshow(Iz);

do3Dplot(Ix,Iy,Iz);
view([90 80 40]);
light('Position',[5 8 12],'Style','local')
axis off
end

