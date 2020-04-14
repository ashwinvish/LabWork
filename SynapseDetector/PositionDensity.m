function  [img] = PositionDensity(cellID,scale,jitter)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

LoadDataFrame;

imageFileName = 'aLZBPointsTm90-37-B45-60-norm-COM-RightPos.tif';
if ismac
    imageFilePath = '/Users/ashwin/Google Drive/Zfish/RefBrains/';
else
    imageFilePath = '/usr/people/ashwinv/seungmount/research/Ashwin/Z-Brain/RefBrains/';
end

imagesInfo = imfinfo(fullfile(imageFilePath,imageFileName));
width = imagesInfo(1).Width;
height = imagesInfo(1).Height;
vol = zeros(height, width,length(imagesInfo));

TifLink = Tiff(fullfile(imageFilePath,imageFileName),'r');
for i = 1:length(imagesInfo)
    TifLink.setDirectory(i);
    vol(:,:,i) = TifLink.read();
end
TifLink.close();

%warning('off','last');
imagesRef = imref3d(size(vol),0.798,0.798,2); % set the correct dimensions for the images using the ref object

Origins = getOrigin(cellID); % origin is in the transformed space

halfX = 2*scale;
halfY = 2*scale;
halfZ = scale;

%img = zeros(halfY*2 + 1,halfX*2 + 1,length(cellID));

randInt = [-1,1];

for jj = 1:length(cellID)
    if sum(Origins(jj,:),2)>0
        cellCenterX = round(Origins(jj,1)/imagesRef.PixelExtentInWorldX)+jitter*randInt(randperm(2,1));
        cellCenterY = round(Origins(jj,2)/imagesRef.PixelExtentInWorldY)+jitter*randInt(randperm(2,1));
        cellCenterZ = round(Origins(jj,3)/imagesRef.PixelExtentInWorldZ)+jitter*randInt(randperm(2,1));
        %[cellCenterX,cellCenterY,cellCenterZ]
        volNans = vol(cellCenterY-halfY:cellCenterY+halfY,cellCenterX-halfX:cellCenterX+halfX,cellCenterZ-halfZ:cellCenterZ+halfZ);
        %volNans(volNans == 0) = NaN;
        img{jj} = volNans;
        %img(:,:,jj) = nanmax(volNans,[],3);
        clear volNans;
    end
end

end

