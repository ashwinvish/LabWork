function parentNode = Parent( tree, node )
%PARENT locates the list of parental nodes
%   tree is the tree structure
%   node is the node whose parents are to be determined

parentNode = tree{node}{1};
for temp = 1:numel(parentNode)
    parentNode = [parentNode, Parent(tree,tree{node}{1}(1))];
end
end

