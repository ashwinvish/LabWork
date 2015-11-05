function [ B, T, BO ] = TreeBranches( tree )
%   TreeBranch computes the number of branches of a tree and the termination
%   points of the same tree

adjMat = AdjMat(tree);                      % Adjacency matrix
typeN = (ones(1,(size(adjMat,1)))*adjMat)';
B = (typeN>1);                             % Branch nodes
T = (typeN==0);                            % Treminal nodes
typeN(typeN>2) = 2;

N = size(adjMat, 1);                        % number of nodes in tree
%typeN = (ones(1,size(dA,1))*dA)';
sdA = adjMat*spdiags(typeN, 0, N, N);
BO = sdA(:,1);                              % calculating weighted path length
resBO = BO;
while sum(resBO)~=0,
    resBO = sdA*resBO;                      % use adjacency matrix to walk through tree
    BO = BO + resBO;
end
BO(1)=1;
BO = full(log2(BO));
end

