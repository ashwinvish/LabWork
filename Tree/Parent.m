function parentNode = Parent_tree( tree, node )
%PARENT locates the list of parental nodes
%   tree is the tree structure
%   node is the node whose parents are to be determined

parentNode = tree{node}{1};
for temp = 1:numel(parentNode)
    parentNode = [parentNode, Parent_tree(tree,tree{node}{1}(1))];
end
end

