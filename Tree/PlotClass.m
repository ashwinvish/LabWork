function [ h ] = PlotClass( cellClass, allTrees, cellIDs, col )
%PLOTCLASS Plots all the cells in a class
%   cellClass is the IDs of all cells in a class, e.g. cellIDsAlx
%   allTrees is the tree structure of all trees
%   cellIDs is the IDs of all cells
%   col is the color of the class

load CellAxons.mat
load LoadSynapses.mat
colors = distinguishable_colors(numel(cellClass),col);
index= 1;

for i = 1:numel(cellIDs)
    if ismember(cellIDs{i}, cellClass) == 1 
        %DisplayTree(allTrees{i},[1],false, [eval([cellIDs{i},'_axon'])],colors(index,:), allPreSynapse{i}, allPost{i});
        DisplayTree(allTrees{i},[1],false, [eval([cellIDs{i},'_axon'])],colors(index,:));
        index = index+1;
    else
        continue
    end
end

set(gca,'Color',[col, 0.2]);
axis vis3d;

end
