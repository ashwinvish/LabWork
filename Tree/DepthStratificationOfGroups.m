function [] = DepthStratificationOfGroups(group, cellIDs, allTrees, col)
%DEPTHSTRATIFICATIONOFGROUPS plots the distribution of synapses along the
%XZ and YZ axes of a group of trees
%   group is cell of strigs with names of trees e.g.cellIDsAlx
%   cellIDS is a cell of strings with the IDs of all the cells
%   allTrees is a cell with the tree structure of all cells in the data
%   col is the 1x3 vector indicating color

% load variables
load CellAxons.mat
load LoadSynapses.mat
groupName = inputname(1);
fname = sprintf('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/ZfishpaperFigures/SupFigures/SupFig6/%sArborDensity.svg',groupName);

% define the varibles
postSynapsesXZ = [];
preSynapsesXZ = [];
postSynapsesYZ = [];
preSynapsesYZ = [];
postSynapsesXY = [];
preSynapsesXY = [];
denXZ = [];
denYZ = [];
denXY = [];


% plot all cells in group
for i = 1:numel(cellIDs)
    if ismember (cellIDs{i}, group)==1
        %h1 = figure(1);
        h1 = subplot(2,2,1);
        DisplayTree(allTrees{i},[1],false, [eval([cellIDs{i},'_axon'])],col,allPreSynapse{i}, allPost{i});
        denTree = [];
        DenTree = allTrees{i};                                                  % iterating through all trees
        validNodes =  eval([cellIDs{i},'_axon']);                   % keeping track of axonal nodes
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
                end
                denTree = [denTree; DnTempx DnTempy DnTempz];            % populate with all dendritic trees of all trees
            end
        end
        
    else
        continue;
    end
    
    edge1 = 0:5:60; % bins of 5um each along the Z axis
    [x1,y1] = histcounts(allPost{i}(:,3)/1000,edge1);
    postSynapsesXZ =  [postSynapsesXZ; [0,x1]];
    [x2,y2] = histcounts(allPreSynapse{i}(:,3)/1000,edge1);
    preSynapsesXZ = [preSynapsesXZ; [0,x2]];
    [x3,y3] = histcounts(-1*denTree(:,3)/1000,edge1);
    denXZ =  [denXZ; [0,x3]];
    
    edge2 = 60:5:240; % bins of 5um each along the Y axis
    [x4,y4] = histcounts(allPost{i}(:,2)/1000,edge2);
    postSynapsesYZ =  [postSynapsesYZ; [0,x4]];
    [x5,y5] = histcounts(allPreSynapse{i}(:,2)/1000, edge2);
    preSynapsesYZ = [preSynapsesYZ; [0,x5]];
    [x6,y6] = histcounts(denTree(:,2)/1000,edge2);
    denYZ =  [denYZ; [0,x6]];
    
    
    edge3 = 20:5:140; % bins of 5um each along the X axis
    [x7,y7] = histcounts(allPost{i}(:,1)/1000,edge3);
    postSynapsesXY =  [postSynapsesXY; [0,x7]];
    [x8,y8] = histcounts(allPreSynapse{i}(:,1)/1000, edge3);
    preSynapsesXY = [preSynapsesXY; [0,x8]];
    [x9,y9] = histcounts(denTree(:,1)/1000,edge3);
    denXY =  [denXY; [0,x9]];
end

% XZ histogram
h2 = subplot(2,2,2);
plot(y1,sum(postSynapsesXZ,1), '-','Color',[0.9,0,0], 'LineWidth',2);
hold on;
[ax,hsub1,hsub2] = plotyy(y1,sum(preSynapsesXZ,1),y3, nanmean(denXZ,1)./ max(nanmean(denXZ,1)))
set(hsub1,'Color',[0,0.8,0], 'LineWidth',2);
set(ax(1),'YColor','k','XDir','reverse','XTick',[0 20 40 60],'YLim',[0, max(sum(postSynapsesXZ,1))],'XaxisLocation','top','FontName', 'Arial', 'FontSize', 40, 'LineWidth',2,'Position',[0.4,0.584,0.335,0.341]);
pbaspect(ax(1),[0.38,0.38,0.38])
set(hsub2,'Color',col, 'LineWidth',2);
set(ax(2),'YColor',col,'YLim',[0 1],'XDir','reverse','YTick',[0.5,1],'XaxisLocation','top','FontName', 'Arial', 'FontSize', 40, 'LineWidth',2,'Position',[0.4,0.584,0.335,0.341]);
view(h2,-90,90);
pbaspect(ax(2),[0.38,0.38,0.38])
box(h2, 'off');
set(h2,'color','none');


%h4 = figure(4);
h4= subplot(2,2,3);
plot(y7 ,sum(postSynapsesXY,1),'-','Color', [0.9,0,0], 'LineWidth',2);
hold on;

[ax,hsub3,hsub4] = plotyy(y8 , sum(preSynapsesXY,1),y9,nanmean(denXY,1)./max(nanmean(denXY,1)))
set(hsub3,'Color',[0,0.8,0], 'LineWidth',2);
set(ax(1),'XLim',[20 140],'YColor','k','YLim',[0 max(sum(postSynapsesXY,1))],'XDir','reverse','XTickLabel',[0,20,40,60,80,100,120],'FontName', 'Arial', 'FontSize', 40, 'LineWidth',4, 'Position',[0.13,0.22,0.335,0.341]);
set(hsub4,'Color',col, 'LineWidth',2);
set(ax(2),'XLim',[20 140],'XDir','reverse','YColor',col,'FontName', 'Arial', 'FontSize', 40, 'LineWidth',4, 'Position',[0.13,0.22,0.335,0.341]);
% relabel the x dim to represent midline as 0; map 140 -->0
box(h4, 'off');
set(h4,'color','none');
% insert stripe locations, 0.293*CellXAxis(loc)
plot(ax(1), [42.19+20, 42.19+20], [0, 200], 'k--', 'LineWidth',4); 
plot(ax(1), [55.08+20, 55.08+20], [0, 200], 'k--','LineWidth',4); 
plot(ax(1), [105.48+20, 105.48+20], [0, 200], 'k--','LineWidth',4); 


% make pretty
view(h1,-180,0);
axis(h1, 'normal');
set(h1,'Color','none','XColor','none','YColor','none', 'ZColor','none');
box(h1,'off');
set(gcf, 'units','normalized','outerposition',[0 0 1 1],'Renderer','painters');



