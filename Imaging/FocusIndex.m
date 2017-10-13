%%

% stdev that straddle the pixel resoultion

imageInfo = imfinfo('TestImage1.tif');
N = imageInfo.Width;
load currentParamaters.mat
load workingDistanceRange.mat

std1 = 6; % upper limit, gaussian with stdev in pixels
std2 = 3; % lower limit, gaussian with stdev in pixels

alpha1 = (N-1)/(2*std1);
alpha2 = (N-1)/(2*std2);


for i = 1:1:21
filename = sprintf('TestImage%d.tif',i);
im1 =  imread(filename);
A = imgaussfilt(im1,alpha1);
B = imgaussfilt(im1,alpha2);
rootMean(i) =  sum(rms(A-B))/size(A,1);
clear im1;
clear A;
clear B;
end

figure; plot(1:1:21,rootMean/max(rootMean),'o');
ylabel('Normalized Focus Index');
[m,n] = max(rootMean);

figure; plot(workingRange*1000,'.');
hold on; plot(n, currentWorkingDistance*1000,'ro');
plot(11, workingRange(11)*1000,'go');
ylabel('Working Distance in mm');