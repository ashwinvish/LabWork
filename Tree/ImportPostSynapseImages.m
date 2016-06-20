% load Tiff stack of all postsynaptic sites

FileTif = '/usr/people/ashwinv/Dropbox/Scripts/PostSynapseROIs.tif';
InfoImage = imfinfo(FileTif);
mImage = InfoImage(1).Width;
nImage = InfoImage(1).Height;
NumberOfImages = length(InfoImage);

PostSynapseImages = zeros(nImage, mImage, NumberOfImages, 'uint8');
TifLink = Tiff(FileTif,'r');

for i = 1:NumberOfImages
    TifLink.setDirectory(i);
    PostSynapseImages(:,:,i) = TifLink.read();
end
TifLink.close();

%%

for i = 1:100
tight_subplot(10,10,0.1,i);
imagesc(PostSynapseImages(:,:,i));
hold on; axis off; axis square
colormap gray;
end
