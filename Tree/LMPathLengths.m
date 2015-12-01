function [LMPathlen] = LMPathLengths(swcPath,rootNode, NodeType)
% LMPATHLENGTHS returns the pathlength of all nodes of type NodeType
%   swcPath is the path to the .swc file
%   rootNode is the location of the first node in the tree, usally 1
%   NodeType is the nodes for which pathlength is required, e.g 2 for axon

[tree,SpecialPositions] = treeVisAV(swcPath,rootNode);
index = 1;
NodeCoords = zeros(numel(tree),3);
for i = 1:numel(tree)
    NodeCoords(i,:) = tree{i}{3};
    if ismember(NodeCoords(i,:),SpecialPositions,'rows')
        NodeCoordsAxon(index,:) = tree{i}{3};
        index = index+1;
    end
end


if strcmp(swcPath,'fish3075_118-axons.swc') == 1
    resolution = [0.36,0.36,1];
else
    resolution = [0.36,0.36,2];
end


LMPathlen = findPathLength(swcPath,tree,resolution,NodeCoordsAxon,NodeType);
