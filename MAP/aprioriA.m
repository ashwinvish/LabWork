

function [ pA ] = aprioriA(A , SigmaA, MuA)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

const = 1./(sqrt(2*pi)*SigmaA) ;
denom = 2*SigmaA.*SigmaA;

for i = 1:1:size(A,3)
    pA(:,i) = const(1)* exp(- (A(:,1,i) - MuA(1)).^2 / denom(1) ) ;
    pA(:,i) = pA(:,i).*(const(2)* exp(- (A(:,2,i) - MuA(2)).^2 / denom(2))) ;
    pA(:,i) = pA(:,i).*(const(3)* exp(- (A(:,3,i) - MuA(3)).^2 / denom(3))) ;
end


pA = pA(:,1)*pA(:,2)';
end

