function [ cellMat ] = AdjMat( tree )
%Adjmat is used to compute the adjacency matrix of a tree
%   The definition of adjaceny matrix is according to the trees toolbox by
%   Kuntz et.al.
%   Created by Ashwin Vishwanathan, ashwinv@princeton.edu, 7/14/2015.

cellMat = zeros(numel(tree));

for kk = 1:numel(tree)
    children = tree{kk}{2};
    for jj = 1:numel(children)
        cellMat(children(jj),kk) = 1;
    end
end

end





