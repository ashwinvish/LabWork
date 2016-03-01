function[lengthToNode] = findPathLength_old(fileName,resolution,queryNodes,validNodeTypes)
% Input queryNodes need to be in pixel coordinates ONLY

if nargin < 4
    validNodeTypes = [0 2 3 4];
end

[nodes,edges,radii,nodeType] = newReadSWCfile(fileName,validNodeTypes);

nodes(:,1) = nodes(:,1)*resolution(1);
nodes(:,2) = nodes(:,2)*resolution(2);
nodes(:,3) = nodes(:,3)*resolution(3);

%specialPositions = nodes(ismember(nodeType,validNodeTypes),:);
h = hist(edges(:,2),[1:size(nodes,1)]); 
irreducibleNodes = union(find(h~=1),1); % 1 is the root node
numelNodes = numel(irreducibleNodes);

for i =  1:size(queryNodes,1)
    tempx = nodes(:,1) - queryNodes(i,1);
    tempy = nodes(:,2) - queryNodes(i,2);
    tempz = nodes(:,3) - queryNodes(i,3);
    
    %TempDiff = abs([tempx,tempy,tempz]);
    TempDiff = sqrt(tempx.^2+tempy.^2+tempz.^2);
    [c,I] = min(TempDiff(:));
    rawLength = 0;
    [lia,locb] = ismember(I,irreducibleNodes);
    
    for kk = 2:locb
        tmpParent = edges(find(edges(:,1)==kk),2) % assume that the edges are ordered pairs: (child, parent)
        path = nodes(kk,:);
        while tmpParent > locb
            newTmpParent = edges(find(edges(:,1)==tmpParent),2);
            path = [path; nodes(tmpParent,:)];
            tmpParent = newTmpParent;
        end
        path = [path; nodes(tmpParent,:)];
        rawPathLengths = sqrt(sum(diff(path,1,1).^2,2));
        rawLength = rawLength + sum(rawPathLengths);
    end
    lengthToNode(i,1) = rawLength;
end

