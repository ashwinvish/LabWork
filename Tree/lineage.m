function nodeList =  lineage(tree, seedNode)
% LINEAGE is the path form the seedNode to the last edge of the tree
%   tree is the tree structure .swc file
%   seedNode is the node for which the children are required.

nodeList = tree{seedNode}{2};
for child = 1:numel(tree{seedNode}{2})
    nodeList = [nodeList lineage(tree, tree{seedNode}{2}(child))];
end

   