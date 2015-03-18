function [] = threeView( volume )
%threeView will genetrate the three view of a volume
%   To generate the three views of a volume XY, XZ, YZ. 

%plot XY
tempXY = max(volume,[],3);
h2 = imagesc(0,90,imrotate(tempXY,-90));set(gca,'YDir','normal'); colormap(jet);
hold on;
axis image;
%plot XZ
tempXZ = max(volume,[],2);
for i = 1:size(tempXZ,3)
    B(:,i) = tempXZ(:,:,i);
end
imagesc(0,0,imrotate(B,-90)); colormap(jet);
axis image;
%plot YZ
tempYZ = max(volume,[],1);
for i = 1:size(tempYZ,3)
    C(i,:) = tempYZ(:,:,i);
end
imagesc(160,90,imrotate(C,-90));set(gca,'YDir','normal'); colormap(jet);
axis image


    
end

