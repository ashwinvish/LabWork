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

% define the varibles
postSynapsesXZ = [];
preSynapsesXZ = [];
postSynapsesYZ = [];
preSynapsesYZ = [];


for i = 1:numel(cellIDs)
    if ismember (cellIDs{i}, group)==1
        %h1 = figure(1);
        h1 = subplot(2,2,1);
        DisplayTree(allTrees{i},[1],false, [eval([cellIDs{i},'_axon'])],col,allPreSynapse{i}, allPost{i});
    else
        continue;
    end
    edge1 = 0:5:60; % bins of 5um each along the XZ axis
    [x1,y1] = histcounts(allPost{i}(:,3)/1000,edge1);
    postSynapsesXZ =  [postSynapsesXZ; [0,x1]];
    [x2,y2] = histcounts(allPreSynapse{i}(:,3)/1000,edge1);
    preSynapsesXZ = [preSynapsesXZ; [0,x2]];
    
    edge2 = 60:5:240; % bins of 5um each along the YZ axis
    [x3,y3] = histcounts(allPost{i}(:,2)/1000,edge2);
    postSynapsesYZ =  [postSynapsesYZ; [0,x3]];
    [x4,y4] = histcounts(allPreSynapse{i}(:,2)/1000, edge2);
    preSynapsesYZ = [preSynapsesYZ; [0,x4]];
end

%h2 = figure(2); % XZ histogram
h2 = subplot(2,2,2);
plot(nanmean(postSynapsesXZ,1)./ max(nanmean(postSynapsesXZ,1)),y1, '-','Color', [0.9,0,0], 'LineWidth',2);
hold on;
plot(nanmean(preSynapsesXZ,1)/max(nanmean(preSynapsesXZ,1)) ,y1,'-','Color',[0,0.8,0], 'LineWidth',2);

%h3 = figure(3);
h3= subplot(2,2,3);
plot(y3 , nanmean(postSynapsesYZ,1)./ max(nanmean(postSynapsesYZ,1)), '-','Color', [0.9,0,0], 'LineWidth',2);
hold on;
plot(y4 , nanmean(preSynapsesYZ,1)./max(nanmean(preSynapsesYZ,1)),'-','Color',[0,0.8,0], 'LineWidth',2);

% make pretty
view(h1,-90,0);
axis(h1, 'normal');
% axis(gca(h1),'vis3d');
set(h1,'Color','none','XColor','none','YColor','none', 'ZColor','none');
box(h1,'off');

pbaspect(h2,[0.38,0.38,0.38])
set(h2,'YLim',[0 60],'XTick',[],'YDir','reverse','XDir','reverse','YTick',[0 20 40 60],'YAxisLocation','right','XColor','none','FontName', 'Arial', 'FontSize', 40, 'LineWidth',2,'Position',[0.4,0.584,0.335,0.341]);
box(h2, 'off');
%axis(gca(h2),'vis3d');
set(h2,'color','none');

%pbaspect(gca(h3),[6,3,1])
set(h3,'XLim',[60 250],'XTick',[60  100  140  180  220  250 ],'YTick',[],'YColor','none','FontName', 'Arial', 'FontSize', 40, 'LineWidth',2, 'Position',[0.13,0.22,0.335,0.341]);
box(h3, 'off');
%axis(gca(h3),'vis3d');
set(h3,'color','none');
set(gcf,'color','w');
