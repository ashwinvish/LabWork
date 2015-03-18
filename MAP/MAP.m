function [ p,A ] = MAP( image1, image2, T1, T2 , AccV, PixSize, n)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%fname1 = image1;
%fname2 = image2;

%im1 = imread(fname1);
im1 = im2double(image1);
%im2 = imread(fname2);
im2 = im2double(image2);

na = NAper(AccV, PixSize);

Sigma = mean([std(im1(:)), std(im2(:))]);

a1 = mean([T1(1),T2(1)]);
a2 = mean([T1(2),T2(2)]);
a3 = mean([T1(3),T2(3)]);

bins =100;

for i = 1:1:n
    
    temp1 = eval(['T' num2str(i) '(1)']);
    temp2 = eval(['T' num2str(i) '(2)']);
    temp3 = eval(['T' num2str(i) '(3)']);
    
    A1 = [temp1-0.1*temp1 :((temp1+0.1*temp1) - (temp1-0.1*temp1))/bins : temp1+0.1*temp1]';
    A2 = [temp2-0.1*temp2 :((temp2+0.1*temp2) - (temp2-0.1*temp2))/bins : temp2+0.1*temp2]';
    A3 = [temp3-0.1*temp3 :((temp3+0.1*temp3) - (temp3-0.1*temp3))/bins : temp3+0.1*temp3]';
    
    A(:,:,i)=[A1,A2,A3];
    
end

MuA = [a1, a2, a3]; SigmaA = [0.05 , 0.01 , 0.01];


pA = aprioriA( A, SigmaA, MuA);
[m,n] = size(im1);

mtf = calculateMTF(A ,na, m,n, PixSize);

fI1 = fft2(im1);
fI2 = fft2(im2);

R = numel(im1);
%constant = 1/(4*pi^2*Sigma^4)^R

tmp = sum(sum( - abs([fI2.*mtf(:,:,1) - fI1.*mtf(:,:,2)]).^2 ./ (2*Sigma*[mtf(:,:,1).^2 + mtf(:,:,2).^2]+1e-10)));

P = log(pA) + tmp;

[p,posp] = max(P)

A = A(posp);

end

