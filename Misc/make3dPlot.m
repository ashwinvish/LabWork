function make3dPlot(s, cmap)

chan = h5read('chunk_49601-50624_28097-29120_17513-17640.img.h5', '/main');
Ix = toColormappedChannel(squeeze(s(1024,:,:)), cmap, squeeze(chan(1024,:,:)), 0.1); imshow(Ix);
Iy = toColormappedChannel(squeeze(s(:,1024,:)), cmap, squeeze(chan(:,1024,:)), 0.1); imshow(Iy);
Iz = toColormappedChannel(squeeze(s(:,:,128)), cmap, squeeze(chan(:,:,128)), 0.1); imshow(Iz);

% Ix = imadjust(Ix,[.2 .2 .2; .7 .7 .7],[]);
% Iy = imadjust(Iy,[.2 .2 .2; .7 .7 .7],[]);
% Iz = imadjust(Iz,[.2 .2 .2; .7 .7 .7],[]);

do3Dplot(Ix,Iy,Iz);
view([90 80 40]);
light('Position',[5 8 12],'Style','local')
axis off

end