function [LMaxLength] = LMPathLengths(swcPath,rootNode, NodeType)
% LMPATHLENGTHS returns the pathlength of all nodes of type NodeType
%   swcPath is the path to the .swc file
%   rootNode is the location of the first node in the tree, usally 1
%   NodeType is the nodes for which pathlength is required, e.g 2 for axon

[tree,SpecialPositions] = treeVisAV(swcPath,rootNode); % SpecialPositions returns all axonal nodes
index = 1;
NodeCoords = zeros(numel(tree),3);
for i = 1:numel(tree)
    NodeCoords(i,:) = tree{i}{3};
    if ismember(NodeCoords(i,:),SpecialPositions,'rows')
        NodeCoordsAxon(index,:) = tree{i}{3};
        AxonalTreeNodes(index) = i;
        index = index+1;
    end
end
LMaxLength = [];
tempLength = 0;
for jj = 1:numel(AxonalTreeNodes)
    tempLength = tempLength + sum(tree{AxonalTreeNodes(jj)}{1,4}{1,2});
end
LMaxLength = [LMaxLength,tempLength];

