function A = SNR(H)
%SNR() Used to estimate the Signal to Noise ratio from single image
%   The noise free ACF is estimated based on JTL Thong et.al. Scanning 2001.

%im = imread(H);
%im = H;
im = im2double(im);
[n,m] = size(im);

im = im/max(max(im)); % normalize image helps to compare


B = autocorr2d(im); % calculate ACF of the image

mu = mean(im(:));
%v = var(im);


B_nf = (B(n/2 + 2,m/2 + 1) + B(n/2 + 1, m/2 + 2)) / 2; % Estimate of ACF of Noise Free image is the average of neighbors

%B_nf = (B(n/2 + 2,m/2 + 1) + B(n/2 + 1, m/2 + 2) + B(n/2 , m/2+1) + B(n/2 + 1, m/2)) / 4;


B_sig = B(n/2 + 1,m/2 + 1); % ACF of signal

Sig = (B_nf - mu.^2) ;

A =  Sig / (B_sig - B_nf);

%A = abs(A);

end

