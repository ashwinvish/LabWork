[tree3003,SpecialPositions3003] = treeVis('fish3003_ch1-exported-000.swc',1);
lin3003 = lineage(tree3003,1,1);
index = 1;
for i = 1:numel(tree3003)
    NodeCoords(i,:) = tree3003{i}{3};
    if ismember(NodeCoords(i,:),SpecialPositions3003,'rows')
        NodeCoordsAxon(index,:) = tree3003{i}{3};
        index = index+1;
    end
end

    length3003 = findPathLength('fish3003_ch1-exported-000.swc',[0.36,0.36,2],NodeCoordAxon,2);
