function sig2noise = SNR(H)
%SNR() Used to estimate the Signal to Noise ratio from single image
%   The noise free ACF is estimated based on JTL Thong et.al. Scanning 2001.
im = H; % read image
im = im2double(im);
[n,m] = size(im); % get width and Height of image

im = im/max(max(im)); % normalize image helps to compare
B = autocorr2d(im); % calculate ACF of the image
mu = mean(im(:)); % mean of image
B_nf = (B(n/2 + 2,m/2 + 1) + B(n/2 + 1, m/2 + 2)) / 2; % Estimate of ACF of Noise Free image as the average of neighbors of center pixel
B_sig = B(n/2 + 1,m/2 + 1); % ACF of signal

Sig = (B_nf - mu.^2) ;
sig2noise =  Sig / (B_sig - B_nf);
end
