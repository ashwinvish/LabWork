function [lengthToNode, SynDiff] = findPathLength_new(fileName,tree,resolution,queryNodes,validNodeTypes)
%FINDPATHLENGTH is to determine the pathlength of the queriedNodes
%   fileName is the .swc file for the tree
%   tree is the tree structre associated with the fileName
%   resolution is the resolution of each pixel [5,5,45]
%   queryNodes are the cartesian coordinates for the nodes that are being
%   queried, need to be in actual resolution.
%   validNodeTypes legacy param.

if nargin < 5
    validNodeTypes = [0 2 3 4];
end


[nodes,edges,radii,nodeType] = newReadSWCfile(fileName,validNodeTypes);

nodes(:,1) = nodes(:,1)*resolution(1);
nodes(:,2) = nodes(:,2)*resolution(2);
nodes(:,3) = nodes(:,3)*resolution(3);


h = hist(edges(:,2),[1:size(nodes,1)]);
irreducibleNodes = union(find(h~=1),1);              % 1 is the root node, irreducible nodes are branch nodes

for kk = 1:numel(irreducibleNodes)
    % swap on edges
    edges(edges==kk) = 0; edges(edges==irreducibleNodes(kk)) = kk; edges(edges==0) = irreducibleNodes(kk);
    % swap on nodes
    tmp = nodes(kk,:); nodes(kk,:) = nodes(irreducibleNodes(kk),:); nodes(irreducibleNodes(kk),:) = tmp;
    
end


for i =  1:size(queryNodes,1)
    tempx = nodes(:,1) - queryNodes(i,1);
    tempy = nodes(:,2) - queryNodes(i,2);
    tempz = nodes(:,3) - queryNodes(i,3);
    
    TempDiff = sqrt(tempx.^2+tempy.^2+tempz.^2);
    [c,I] = min(TempDiff(:));                       % find the queryNode among nodes
    rawLength = 0;
    [m,n] = min(abs(I-irreducibleNodes));           % find closest irreducibleNode
    
    lengthToClosesIrreducibleNode = 0;
    tmpParent = edges(find(edges(:,1)==I),2);
    pth = nodes(I,:);
    
    
    while tmpParent ~= 1                            % parse through all edges, from edge with queryNode to root node 
        newTmpParent = edges(find(edges(:,1)==tmpParent),2);
        pth = [pth; nodes(tmpParent,:)];
        tmpParent = newTmpParent;
    end
    
    pth = [pth; nodes(I,:)];                        % all linesegments till query node
    rawPathLengths = sqrt(sum(diff(pth,1,1).^2,2));
    lengthToNode(i,1) = lengthToClosesIrreducibleNode + sum(rawPathLengths);
    
end

% inter-synaptic distance
SynDiff = diff(sort(lengthToNode));

end




