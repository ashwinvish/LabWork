clear all;
cd('/usr/people/ashwinv/seungmount/research/Ashwin/Zfish/SectionOverviewsAlignedWithTemplateDirectory');
saveDir = '/usr/people/ashwinv/seungmount/research/Ashwin/Zfish/Cbllm';
list = dir('*.tif');
numFiles =  size(list,1);

%preallocate
im1 = zeros(4096,4096);
im2 = zeros(4096,4096);

tic;

for i = 2:numFiles
    im1 = im2double(imread(list(i-1).name));
    im2 = im2double(imread(list(i).name));
    
    [out,greg] = dftregistration(fft2(im1),fft2(im2),100);
    OutFileName = sprintf('AlignedSection_%d.tif',i);
    imwrite(ifft2(greg),fullfile(saveDir,OutFileName));
    
end

toc
