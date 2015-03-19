function [] = threeView( volume, cmap , trans)
%threeView will genetrate the three view of a volume
%   To generate the three views of a volume XY, XZ, YZ. 

%plot XY
tempXY = max(volume,[],3);
h1 = imagesc(0,90,imrotate(tempXY,-90));
set(h1,'AlphaData', trans);
set(gca,'YDir','normal'); 
colormap(cmap);
hold on;
axis image;

%plot XZ
tempXZ = max(volume,[],2);
for i = 1:size(tempXZ,3)
    B(:,i) = tempXZ(:,:,i);
end
h2 = imagesc(0,0,imrotate(B,-90));
set(h2,'AlphaData',trans);
colormap(cmap);
axis image;

%plot YZ
tempYZ = max(volume,[],1);
for i = 1:size(tempYZ,3)
    C(i,:) = tempYZ(:,:,i);
end
h3 = imagesc(160,90,imrotate(C,-90));
set(h3,'AlphaData',trans);
set(gca,'YDir','normal'); 
colormap(cmap);
axis image

axis off;
axis vis3d;
end

