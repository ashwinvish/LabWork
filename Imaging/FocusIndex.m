%%

% stdev that straddle the pixel resoultion

imageInfo = imfinfo('TestImage1.tif');
N = imageInfo.Width;
load currentParamaters.mat
load workingDistanceRange.mat
numImages = 17;

std1 = 6; % upper limit, gaussian with stdev in pixels
std2 = 3; % lower limit, gaussian with stdev in pixels

alpha1 = (N-1)/(2*std1);
alpha2 = (N-1)/(2*std2);


for i = 1:numImages
filename = sprintf('TestImage%d.tif',i);
im1 =  imread(filename);
A = imgaussfilt(im1,alpha1);
B = imgaussfilt(im1,alpha2);
rootMean(i) =  sum(rms(A-B))/size(A,1);
clear im1;
clear A;
clear B;
end


pol = polyfit(workingRange,rootMean,3); % atleast 5 variables, Energy, DwellTime, Beam current, Stigmations x,y

x1= linspace(min(workingRange)*1000,max(workingRange)*1000,10);
y1 = polyval(pol,x1);

figure();
clear p;

p(1,1) =  gramm('x',[workingRange*1000,x1], 'y',[rootMean,y1], 'color',[ones(size(workingRange,2),1)',2*ones(size(x1,2),1)']);
p(1,1).geom_point();
p(1,1).geom_line();
p(1,1).set_names('x','WorkingDistance','y','Focus Index');


[m,n] = max(rootMean);

p(1,2) = gramm('x',[1:numImages,n,9]', 'y', [workingRange*1000,currentWorkingDistance*1000,workingRange(9)*1000]', ...
    'color',[ones(size(1:numImages,2),1);2;3]);
p(1,2).geom_point();
p(1,2).set_names('x','Image Number','y','Working distance');

p.draw();
