function list =  lineage(tree, seedNode, list)

list = [list seedNode];
children = tree{seedNode}{2};
if ~isempty(children)
    for kk=1:numel(children)
      list = lineage(tree,children(kk),list);
    end
end


   