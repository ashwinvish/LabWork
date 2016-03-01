function [lengthToNode, SynDiff] = findPathLength_new(fileName,tree,resolution,queryNodes,validNodeTypes)
%FINDPATHLENGTH is to determine the pathlength of the queriedNodes
%   fileName is the .swc file for the tree
%   tree is the tree structre associated with the fileName
%   resolution is the resolution of each pixel [5,5,45]
%   queryNodes are the cartesian coordinated for the nodes that are being
%   queried, need to be in actual resolution.
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
    
    TempDiff = sqrt(tempx.^2+tempy.^2+tempz.^2);
    [c,I] = min(TempDiff(:));
    rawLength = 0;
    [m,n] = min(abs(I-irreducibleNodes));
    %display(n);
    parentNodes = Parent(tree,n);
    parentNodes = [n,parentNodes];
    PntNodes(i,1:size(parentNodes,2)) = parentNodes;
    for kk = 2:I
        tmpParent = edges(find(edges(:,1)==kk),2); % assume that the edges are ordered pairs: (child, parent)
        path = nodes(kk,:);
        path = [path; nodes(tmpParent,:)];
        rawPathLengths = sqrt(sum(diff(path,1,1).^2,2));
        rawLength = rawLength + sum(rawPathLengths);
    end
    lengthToNode(i,1) = rawLength;
end 

% inter synaptic distance for synapse that share the same parent
    for i = 2:size(queryNodes,1)
        if PntNodes(i-1,1) == PntNodes(i,1)
            SynDiff(i-1) = abs(lengthToNode(i)- lengthToNode(i-1));
        else
            SynDiff(i-1) = 0;
        end
        
    end
    



