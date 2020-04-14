% load images
 imageFilePath = '/Users/ashwin/Google Drive/Zfish/RefBrains/';
imageFileName = 'Alex-v2-RightwardPosition-Blacks.tif';
imagesInfo = imfinfo(fullfile(imageFilePath,imageFileName));
width = imagesInfo(1).Width;
height = imagesInfo(1).Height;
%images = zeros(height, width,length(imagesInfo),'uint8');
imagesPos = zeros(height, width,3,length(imagesInfo));
 TifLink = Tiff(fullfile(imageFilePath,imageFileName),'r');    
    for i = 1:length(imagesInfo)
        TifLink.setDirectory(i);
        imagesPos(:,:,:,i) = TifLink.read();
    end
  TifLink.close();
  
  
  
 imageFileName = 'Alex-v2-RightwardVelocity-Reds.tif';
imagesInfo = imfinfo(fullfile(imageFilePath,imageFileName));
width = imagesInfo(1).Width;
height = imagesInfo(1).Height;
%images = zeros(height, width,length(imagesInfo),'uint8');
imagesVel = zeros(height, width,3,length(imagesInfo));
  TifLink = Tiff(fullfile(imageFilePath,imageFileName),'r');   
    for i = 1:length(imagesInfo)
        TifLink.setDirectory(i);
        imagesVel(:,:,:,i) = TifLink.read();
    end
  TifLink.close();
  
% equlaize histograms

imagesVelmatch = imhistmatchn(imagesVel,imagesPos);
  
  % P-V/P+V
  

  
A = imagesPos;
B = imagesVelmatch;
C = imdivide(A,B);

  
subplot(1,3,1)
imagesc(A(:,:,50))
subplot(1,3,2)
imagesc(B(:,:,50))
subplot(1,3,3)
imagesc(C(:,:,50))