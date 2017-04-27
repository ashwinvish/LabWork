function [ combImage ] = combinedImage( originalImage, alteredImage, downsample )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

if nargin<3
    downsample = 1;
end

originalImage = imresize(originalImage, downsample);
alteredImage = imresize(alteredImage, downsample);

[height,width] =  size(originalImage);
%preallocate
combImage = zeros(width,width,3);

combImage(:,:,1) = originalImage; %red
combImage(:,:,2) = alteredImage; %green

imshow(combImage);

end

