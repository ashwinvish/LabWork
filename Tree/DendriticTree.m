function [denTree] = DendriticTree(tree,treeno, cellIDs, Display)
% DENDRITICTREE is to extract all the nodes of the dendrites of a tree
%   Tree is the tree whos dendrites are being extracted
%   treeno is the ID associated with the tree
%   cellIDs are the IDs of all cells in the dataset
%   Display true or false


load CellAxons.mat

denTree = [];
DenTree = tree;                                                 % iterating through all trees
validNodes =  eval([cellIDs{treeno},'_axon']);                   % keeping track of axonal nodes
if isempty(validNodes)                                           % if no axon check
    validNodes = 1:numel(DenTree);
end                                                              % if no axon then iterate through all nodes
for jj = 1:numel(DenTree)
    children = DenTree{jj}{2};                                   % Consider childern of jj node
    for nn = 1:numel(children)                                   % iterate over all children
        DnTempx = [DenTree{children(nn)}{3}(1); DenTree{children(nn)}{4}{1}(:,1)];
        DnTempy = [DenTree{children(nn)}{3}(2); DenTree{children(nn)}{4}{1}(:,2)];
        DnTempz = -[DenTree{children(nn)}{3}(3); DenTree{children(nn)}{4}{1}(:,3)];
        if length(validNodes)< length(DenTree) && ismember(jj,validNodes) % dont plot the axon
            continue;
        else
            if Display == true
                h2 = plot3(DnTempx,DnTempy,DnTempz,'-','color',[0.8 0.8 0.8]);
                %h2.Color(4) = 0.2;
            end
        end
        denTree = [denTree; DnTempx DnTempy DnTempz];            % populate with all dendritic trees of all trees
    end
end
end
