% AxonCollateral lengths for trees

% Int1_5
Int1_5col{1} = [155,154]
Int1_5col{2} = [152,155]
Int1_5col{3} = [152,156,157,158]
Int1_5col{4} = [140,143]
Int1_5col{5} = [132,134]
Int1_5col{6} = [128,131,136,141,142,144]
Int1_5col{7} = [103,117,123,125,126,127,137,139,150,151,153]
Int1_5col{8} = [124,117,131,132,140,152]
Int1_5col{9} = [92,99,100,103,115]
Int1_5col{10} = [77,92,106,110,114,116,133,135,145,146,147,148,149]
Int1_5col{11} = [77,80,81,82,84,85,86]
Int1_5col{12} = [163,165,168]
Int1_5col{13} = [163,165,166,167,169,170]

LR_Caudal = [56,107,159];
LR_RostrAL = [56,77];

Nodes = [];
for i = 1:numel(Int1_5col)
    axonNodes = AxonQueryNodes(allTrees{5},Int1_5col{i});
    tmpNodes =[];
    for jj = 1:numel(Int1_5col{i})
        tmpParent = Parent(allTrees{5},Int1_5col{i}(jj));
        for kk = 1:numel(tmpParent)
            tmpNodes = [tmpNodes;allTrees{5}{tmpParent(kk)}{1,4}{1,1}];
        end
    end
    Nodes = [Nodes;unique(tmpNodes,'rows')];
    Nodes = unique(Nodes,'rows');
    colPathLengthTemp{i} = findPathLength('Int1_5_WithTags.swc',allTrees{5},[5,5,45],axonNodes);
    Int1_5colPathLength{i} = sum(diff(sort(colPathLengthTemp{i})));
end
Int1_5colPreSynapses = sum(ismember(Nodes,allPreSynapse{5},'rows'));
clear axonNodes;
clear colPathLengthTemp;