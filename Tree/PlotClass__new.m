function [ h ] = PlotClass_new( cellClass, allTrees, cellIDs,rho )
%PLOTCLASS Plots all the cells in a class
%   cellClass is the IDs of all cells in a class, e.g. cellIDsAlx
%   allTrees is the tree structure of all trees
%   cellIDs is the IDs of all cells
%   col is the color of the class

load CellAxons.mat
load LoadSynapses.mat
col = parula(22);
[y,I] = sort(rho);


for i = 1:numel(cellIDs)
    if ismember(cellIDs{i}, cellClass) == 1 
        Pos = find(I==i);
        %DisplayTree(allTrees{i},[1],false, [eval([cellIDs{i},'_axon'])],col(Pos,:), allPreSynapse{i}, allPost{i});
        DisplayTree(allTrees{i},[1],false, [eval([cellIDs{i},'_axon'])],col(Pos,:));

    else
        continue
    end
end

end

