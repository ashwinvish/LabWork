function [ B, T ] = TreeBranches( tree )
%TreeBranch computes the number of branches of a tree and the termination
%   points of the same tree

adjMat = AdjMat(tree);
temp = ones(1,(size(adjMat,1)))*adjMat;
B = (temp>1)';
T = (temp==0)';

end

