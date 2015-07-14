% plot all trees without displaying annotated features
% created on 06/09/2015 

D = dir('*.swc'); % choose all .swc files in the current directory
for i = 1:length(D)
    [path,name,ext] = fileparts(D(i).name);
    fname = name(1:11);
    [fname,rawlength,temp] = generateIrreducibleDoubleLinkedTree_WithDim([name,ext],[-1:10],6,true);
    allTrees{i} = fname; allNames{i} = name;
    treeVisualizer(allTrees{i},[1],[],[],false,{[rand rand rand]}, 1:numel(allTrees{i}), false)
end

