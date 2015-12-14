% test for tree 1_4

AxnNodes1_4 = eval([cellIDs{4},'_axon']);
TerminalNodes1_4 = find(Terminals{4}==0);
CommonNodes1_4 = AxnNodes1_4(find(ismember(AxnNodes1_4,TerminalNodes1_4)==1));
for jj = 1:numel(CommonNodes1_4)
    treeNodes1_4(jj,1:3) = allTrees{4}{CommonNodes1_4(jj)}{1,3};
end
allLengthToPreNode1_4 = findPathLengthNew([cellIDs{4} , '_WithTags.swc'],allTrees{4},[5,5,45],treeNodes1_4);





Int1_4col{1} = [67,68,70]
Int1_4col{2} = [62,65,66,68,71,72,73,74]
Int1_4col{3} = [62,63,65,69]
Int1_4col{4} = [48,49,56]
Int1_4col{5} = [48,51]

Int1_4LR_Rostral = [41,48]
Int1_4LR_Caudal = [41,51,52,53,54,55,57,58,59,60,61]

for i = 1:numel(Int1_4col)
    axonNodes = AxonQueryNodes(allTrees{4},Int1_4col{i});
    colPathLengthTemp{i} = findPathLengthNew('Int1_4_WithTags.swc',allTrees{4},[5,5,45],axonNodes);
    Int1_4colPathLength{i} = sum(diff(sort(colPathLengthTemp{i})));
end

clear axonNodes;
clear colPathLengthTemp;