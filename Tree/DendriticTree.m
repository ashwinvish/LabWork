function [denTree,x,y] = DendriticTree(tree,treeno, cellIDs, col, Display)
% DENDRITICTREE is to extract all the nodes of the dendrites of a tree
%   Tree is the tree whos dendrites are being extracted
%   treeno is the ID associated with the tree
%   cellIDs are the IDs of all cells in the dataset
%   col is a string assocaited with the color. eg. 'b'
%   Display true or false


load CellAxons.mat
denTree = [];
DenTree = tree;                                                  % iterating through all trees
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
                subplot(1,2,1)
                h2 = plot3(DnTempx,DnTempy,DnTempz,'-','color',col);
                hold on;
                %h2.Color(4) = 0.2;
            end
        end
        denTree = [denTree; DnTempx DnTempy DnTempz];            % populate with all dendritic trees of all trees
    end
end
scatter3(denTree(1,1), denTree(1,2),denTree(1,3),500,'MarkerFaceColor',col,'MarkerEdgeColor','k');
daspect([1,1,1]);
view(-180,0);
set (gca,'XTick',[], 'YTick',[],'ZTick', [], 'ZLim',[-60000, -0],'XLim',[20000 , 140000]);
box on;

subplot(1,2,2)
h = histogram(-1*denTree(:,3)/1000, 'BinWidth',1,'Orientation','horizontal','FaceColor',col);
 x= h.Values;
 y = h.BinEdges;
% plot([0,x],y,'-','Color',col);
pbaspect([2,1,1]);
%plot(h.Values,1:2:60-1);
hold on;
set(gca,'YLim',[0 ,60],'YDir','reverse', 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
set(gcf,'color','w');
box off;
end
