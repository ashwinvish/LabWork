clc;
clear all;

 StackDir = uigetdir('','Pick the montage stack directory');
 
% DirAndWildCardStr = sprintf('10',StackDir);
figno = 1244;
 
ListOfFiles = dir([StackDir '/*.tif']);
 
 for i = 1:length(ListOfFiles)
     
     fname = fullfile (StackDir , ListOfFiles(i).name);
     imtemp = imread(fname);
%      maxgray = max(max(imtemp));
%      imtemp = imtemp./maxgray;
     disp(sprintf('Loading File: %s',fname));
     SigNoise(i) = SNR(imtemp);
     %Fourier(i) = fft2(imtemp);
%      temp = fft2(imtemp);
     %figure (figno);
     %subplot(4,4,i); imagesc(abs(fftshift(log10(temp))));
%      colormap gray;
%      imtemp = im2double(imtemp);
%      tempprime = imtemp(:)';
%      Norm(i) = norm(tempprime);
     
 end
 %%
%figure;
hold all;
DwellTime = 0.1:0.3:2;
plot(DwellTime , SigNoise , '-o');

 