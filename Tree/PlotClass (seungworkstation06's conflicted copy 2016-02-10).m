function [ h ] = PlotClass( cellClass, allTrees, cellIDs, bgCol )
%PLOTCLASS Plots all the cells in a class
%   cellClass is the IDs of all cells in a class, e.g. cellIDsAlx
%   allTrees is the tree structure of all trees
%   cellIDs is the IDs of all cells
%   col is the color of the class

load CellAxons.mat
load LoadSynapses.mat
col =   distinguishable_colors(numel(cellClass),bgCol);
index = 1;
for i = 1:numel(cellIDs)
    if ismember(cellIDs{i}, cellClass) == 1
        DisplayTree(allTrees{i},[1],false, [eval([cellIDs{i},'_axon'])],col(index,:));
       % DisplayTree(allTrees{i},[1],false, [eval([cellIDs{i},'_axon'])],col, allPreSynapse{i}, allPost{i});
       index = index+1;
    else
        continue
    end
end
set(gca,'Color',[bgCol,0.2]);
end

