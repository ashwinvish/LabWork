function[lengthToNode] = findPathLength(fileName,resolution,synapticNodes,validNodeTypes)
% Input synapticNodes need to be in pixel coordinates ONLY

if nargin < 4
    validNodeTypes = [0 2 3 4];
end

[nodes,edges] = newReadSWCfile(fileName,validNodeTypes);
nodes(:,1) = nodes(:,1)*resolution(1);
nodes(:,2) = nodes(:,2)*resolution(2);
nodes(:,3) = nodes(:,3)*resolution(3);

for i =  1:size(synapticNodes,1)
    tempx = nodes(:,1) - synapticNodes(i,1);
    tempy = nodes(:,2) - synapticNodes(i,2);
    tempz = nodes(:,3) - synapticNodes(i,3);
    
    %TempDiff = abs([tempx,tempy,tempz]);
    TempDiff = sqrt(tempx.^2+tempy.^2+tempz.^2);
    [c,I] = min(TempDiff(:));
    rawLength = 0;
    
    for kk = 2:I
        tmpParent = edges(find(edges(:,1)==kk),2); % assume that the edges are ordered pairs: (child, parent)
        path = nodes(kk,:);
        while tmpParent > mean(I)
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

