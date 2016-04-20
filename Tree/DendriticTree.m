function [denTree, denXY, denYZ, denXZ, edge1, edge2, edge3] = DendriticTree(tree,treeno, cellIDs, col, Display)
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
                h1 = subplot(2,2,1);
                plot3(DnTempx,DnTempy,DnTempz,'-','color',col);
                hold on;
            end
        end
        denTree = [denTree; DnTempx DnTempy DnTempz];            % populate with all dendritic trees of all trees
    end 
end
% plot somata
plot3(denTree(1,1), denTree(1,2),denTree(1,3),'o','MarkerSize',20,'MarkerFaceColor',col,'MarkerEdgeColor','k');
hold on;

denXZ = [];
denYZ = [];
denXY = [];

edge1 = 0:5:60; % bins of 5um each along the Z axis
[x1,y1] = histcounts(-1*denTree(:,3)/1000,edge1);
denXZ =  [denXZ; [0,x1]];

edge2 = 60:5:240; % bins of 5um each along the Y axis
[x2,y2] = histcounts(denTree(:,2)/1000,edge2);
denYZ =  [denYZ; [0,x2]];

edge3 = 20:5:140; % bins of 5um each along the X axis
[x3,y3] = histcounts(denTree(:,1)/1000,edge3);
denXY =  [denXY; [0,x3]];



if Display == true
    view(h1, -180,0);
    axis(h1, 'normal');
    set (h1,'XTick',[], 'YTick',[],'ZTick', [],'XColor','none','YColor','none', 'ZColor','none','ZLim',[-60000, -0],'XLim',[20000 , 140000]);
    box(h1,'off');
    
    h2 = subplot(2,2,2);
    pbaspect(h2,[0.38,0.38,0.38])
    hold on;
    plot(denXZ./ max(denXZ),y1, '-','Color',col, 'LineWidth',2);
    set(h2, 'YLim',[0 60],'XTick',[0.5,1],'YDir','reverse','XDir','reverse','YTick',[0 20 40 60],'YAxisLocation','right','FontName', 'Arial', 'FontSize', 40, 'LineWidth',2,'Position',[0.4,0.584,0.335,0.341]);
    box(h2, 'off');
    set(h2,'color','none');
    
    h3 = subplot(2,2,3);
    hold on;
    plot(y3,denXY./max(denXY),'-','Color',col, 'LineWidth',2);
    set(h3, 'XLim',[20 140],'XDir','reverse','XTickLabel',[0,20,40,60,80,100,120],'YTick',[0.5,1],'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2 , 'Position',[0.13,0.22,0.335,0.341])
    box (h3, 'off');
    set(h3,'Color','none');
end;


end
