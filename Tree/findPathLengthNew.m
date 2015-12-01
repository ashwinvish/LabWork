function [lengthToNode ] = findPathLength(fileName,tree,resolution,queryNodes,validNodeTypes)
%FINDPATHLENGTH is to determine the pathlength of the queriedNodes
%   fileName is the .swc file for the tree
%   tree is the tree structre associated with the fileName
%   resolution is the resolution of each pixel [5,5,45]
%   queryNodes are the cartesian coordinated for the nodes that are being
%   queried
%   validNodeTypes legacy param

if nargin < 5
    validNodeTypes = [0 2 3 4];
end

[nodes,edges,radii,nodeType] = newReadSWCfile(fileName,validNodeTypes);
h = hist(edges(:,2),[1:size(nodes,1)]);
irreducibleNodes = union(find(h~=1),1); % 1 is the root node

nodes(:,1) = nodes(:,1)*resolution(1);
nodes(:,2) = nodes(:,2)*resolution(2);
nodes(:,3) = nodes(:,3)*resolution(3);

for i =  1:size(queryNodes,1)
    tempx = nodes(:,1) - queryNodes(i,1);
    tempy = nodes(:,2) - queryNodes(i,2);
    tempz = nodes(:,3) - queryNodes(i,3);
    
    %TempDiff = abs([tempx,tempy,tempz]);
    TempDiff = sqrt(tempx.^2+tempy.^2+tempz.^2);
    [c,I] = min(TempDiff(:));
    rawLength = 0;
    %[lia,locb] = ismember(I,irreducibleNodes);    % locate closest tree node
    [m,n] = min(abs(I-irreducibleNodes));
    %display(n);
    parentNodes = Parent(tree,n);
    parentNodes = [n,parentNodes];
    
    for j = 1:numel(parentNodes)
        tempParent = parentNodes(j);
        PathLengthToParent = sum(tree{tempParent}{1,4}{1,2});
        rawLength = rawLength+PathLengthToParent;
    end
    lengthToNode(i,1) = rawLength;
end
end

