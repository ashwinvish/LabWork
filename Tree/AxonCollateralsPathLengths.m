% AxonCollateral lengths for trees

% Int1_5
Int1_5col{1} = [117,124,128,131,132,134,136,138,140,141,142,143,144,152,154,155,156,157,158];
Int1_5col{2} = [103,117,123,125,126,127,137,139,150,151,153];
Int1_5col{3} = [92,99,100,103,115];
Int1_5col{4} = [77,92,106,110,114,116,133,135,145,146,147,148,149];
Int1_5col{5} = [159,161,163];
Int1_5col{6} = [163,165,168];
Int1_5col{7} = [163,165,166,167,169,170];
Int1_5col{8} = [77,80,81,82,84,85,86];

LR_Caudal = [56,107,159];
LR_RostrAL = [56,77];


for i = 1:numel(Int1_5col)
    axonNodes = AxonQueryNodes(allTrees{5},Int1_5col{i});
    colPathLengthTemp{i} = findPathLength('Int1_5_WithTags.swc',allTrees{5},[5,5,45],axonNodes);
    Int1_5colPathLength{i} = sum(diff(colPathLengthTemp{i}));
end

%%
% Int1_4
Int1_4col{1} = [67,68,70]
Int1_4col{2} = [62,65,66,68,71,72,73,74]
Int1_4col{3} = [62,63,65,69]
Int1_4col{4} = [48,49,56]
Int1_4col{5} = [48,51]

Int1_4LR_Rostral = [41,48]
Int1_4LR_Caudal = [41,51,52,53,54,55,57,58,59,60,61]

for i = 1:numel(Int1_4col)
    axonNodes = AxonQueryNodes(allTrees{4},Int1_4col{i});
    colPathLengthTemp{i} = findPathLength('Int1_4_WithTags.swc',allTrees{4},[5,5,45],axonNodes);
    Int1_4colPathLength{i} = sum(diff(colPathLengthTemp{i}));
end
 
%%
% Int1_6

Int1_6col{1} = [130,138]
Int1_6col{2} = [130,137,145]
Int1_6col{3} = [137,147]
Int1_6col{4} = [74,118,130]
Int1_6col{5} = [136,146];
Int1_6col{6} = [74,118,120,121,136]
Int1_6col{7} = [120,136,162]
Int1_6col{8} = [162,168]
Int1_6col{9} = [162,169,170,171,172,173]
Int1_6col{10} = [148,155]
Int1_6col{11} = [148,150,151,157]
Int1_6col{12} = [144,149,152,153]
Int1_6col{13} = [144,134,148]
Int1_6col{14} = [134,139]
Int1_6col{15} = [123,126,128,133]
Int1_6col{16} = [131,140,142,156,158,160,161,164,163,159,165,167]
Int1_6col{17} = [131,135,166]
Int1_6col{28} = [135,141,143,154]


Int1_6LR_Caudal = [82,96,98,99,134]
Int1_6LR_Rostral = [82,123,129,131]

for i = 1:numel(Int1_6col)
    axonNodes = AxonQueryNodes(allTrees{6},Int1_6col{i});
    colPathLengthTemp{i} = findPathLength('Int1_6_WithTags.swc',allTrees{6},[5,5,45],axonNodes);
    Int1_6colPathLength{i} = sum(diff(colPathLengthTemp{i}));
end

%%
% Int1_7

Int1_7col{1} = [115,129]
Int1_7col{2} = [115,120]
Int1_7col{3} = [81,115,99]
Int1_7col{4} = [30,75,81]
Int1_7col{5} = [75,110,114]
Int1_7col{6} = [110,116,117,121,122,123]

Int1_7LR_Contra = [30,76]

for i = 1:numel(Int1_7col)
    axonNodes = AxonQueryNodes(allTrees{7},Int1_7col{i});
    colPathLengthTemp{i} = findPathLength('Int1_7_WithTags.swc',allTrees{7},[5,5,45],axonNodes);
    Int1_7colPathLength{i} = sum(diff(colPathLengthTemp{i}));
end

%%






  
