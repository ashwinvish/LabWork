function list =  lineage(tree, seedNode,list)
% Lineage is the path form the seedNode to the last edge of the tree

list = [list seedNode];
children = tree{seedNode}{2};
if ~isempty(children)
    for kk=1:numel(children)
      list = lineage(tree,children(kk),list);
    end
end


   