function [mtf] = calculateMTF( A, na, m, n, PixSize )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
mtf = zeros(m,n,size(A,3));

[Kx,Ky] = meshgrid ([1:n]/(PixSize*n*1e-6),[1:m]/(PixSize*m*1e-6) );



% for i = 1:1:n
%     for j = 1:1:m
% 
%         Kx = i/(PixSize*n*1e-6);
%         Ky = j/(PixSize*m*1e-6);

for k = 1:1:size(A,3)
     %mtf(i,j,k) = exp(-0.125*(na^2)*((2*A(1,2,k)*(Kx^2 - Ky^2)*A(1,1,k) - A(1,3,k).*Kx*Ky*A(1,1,k) + (Kx^2+Ky^2)*(A(1,2,k).^2 + A(1,3,k)^2 + A(1,1,k)^2))));
    
    mtf(:,:,k) = exp(-0.125*(na^2)*((2*A(1,2,k)*(Kx^2 - Ky^2).*A(1,1,k) - A(1,3,k).*Kx*Ky*A(1,1,k) + (Kx^2+Ky^2).*(A(1,2,k).^2 + A(1,3,k)^2 + A(1,1,k)^2))));
end


%     end
% end

