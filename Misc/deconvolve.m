% decolvolve images series
% created on 3/16/2012 by Ashwin

clc;
clear all;

fname = '101112_1_023.tif';
info = imfinfo(fname);
num_images = numel(info);

width = info.Width;
height = info.Height;

J = zeros(height, width, num_images); % allocating space for deconv images
I= zeros(height, width, num_images); % allocating space for input
psf = zeros(height, width, num_images); % allocating space for psfs
Inipsf = ones(height, width);

for i = 1: num_images
    I(:,:,i) = imread(fname, i, 'Info', info);
    [J(:,:,i), psf(:,:,1)] = deconvblind(I(:,:,i), Inipsf);
    J(:,:,i) = J(:,:,i)-min(min(J(:,:,i))); J(:,:,i) = J(:,:,i)/max(max(J(:,:,i)));
    J(:,:,i) = round(J(:,:,i)*255);
    imwrite(uint16(J(:,:,i)), 'decon_Zfish.tif', 'WriteMode', 'append', 'Compression', 'none');
    imwrite(uint16(psf(:,:,i)), 'decon_psf_Zfish.tif', 'WriteMode', 'append', 'Compression', 'none');
    disp(i)
end


