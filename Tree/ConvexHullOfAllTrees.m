TotalVolume = 120000*190000*60000; % in nm


AlxHull = [];
AlxVol = [];
TransHull = [];
TransVol = [];
DbxHull = [];
DbxVol = [];
BarhlHull = [];
BarhlVol =[];


for i = 1:22
    denTree = DendriticTree(allTrees{i},i,cellIDs,calx,false);
        [ConHull ConVol] = convhull(denTree(:,1),denTree(:,2),denTree(:,3),'simplify', true);
        htri1 = trimesh(ConHull,denTree(:,1),denTree(:,2),denTree(:,3),'faceColor',CellColor(i,:),'FaceAlpha',0.2,'EdgeColor','black','EdgeAlpha',0.1);
        hold on;
        scatter3(CellSoma(i,1),CellSoma(i,2),CellSoma(i,3),200, 'Marker','o', 'MarkerFaceColor',CellColor(i,:), 'MarkerEdgeColor','k','LineWidth', 2);
        daspect([1 1 1]);
        clear ConHull;
        clear ConVol;
        clear denTree;
        pause(1);
end

%     axis vis3d;
%         axis([20000 140000  60000 250000 -60000 0]);
%         plot( [20000, 40000], [70000, 70000],'-k' ) % insert 20um sclaebar
%         box off;
%         XColor = [1,1,1]; YColor = [1,1,1];
%         set (gca, 'XTick',[], 'YTick',[],'ZTick', [], 'Ydir','reverse');
%         view([-180,90]);

% 
% for i = 1:numel(cellIDs)
%     if ismember(cellIDs{i}, cellIDsAlx) ==1
%         figure(1);
%         denTree = DendriticTree(allTrees{i},i,cellIDs,calx,false);
%         [ConHull ConVol] = convhull(denTree(:,1),denTree(:,2),denTree(:,3),'simplify', true);
%         AlxHull{i} = ConHull;
%         AlxVol = [AlxVol, ConVol];
%         htri1 = trimesh(ConHull,denTree(:,1),denTree(:,2),denTree(:,3),'faceColor',calx,'FaceAlpha',0.2,'EdgeColor','black','EdgeAlpha',0.1);
%         hold on;
%         treeCol(i,:) = calx;
%         scatter3(CellSoma(i,1),CellSoma(i,2),CellSoma(i,3),200, 'Marker','o', 'MarkerFaceColor',calx, 'MarkerEdgeColor','k','LineWidth', 2);
%         daspect([1 1 1]);
%         axis vis3d;
%         axis([20000 140000  60000 250000 -60000 0]);
%         plot( [20000, 40000], [70000, 70000],'-k' ) % insert 20um sclaebar
%         box off;
%         XColor = [1,1,1]; YColor = [1,1,1];
%         set (gca, 'XTick',[], 'YTick',[],'ZTick', [], 'Ydir','reverse');
%         view([-180,90]);
%         clear ConHull;
%         clear ConVol;
%         clear denTree;
%     elseif ismember(cellIDs{i}, cellIDsTrans) ==1
%         figure(2);
%         denTree = DendriticTree(allTrees{i},i,cellIDs,ctrans,false);
%         [ConHull ConVol] = convhull(denTree(:,1),denTree(:,2),denTree(:,3),'simplify', true);
%         TransHull{i} = ConHull;
%         TransVol = [TransVol, ConVol];
%         htri1 = trimesh(ConHull,denTree(:,1),denTree(:,2),denTree(:,3),'faceColor',ctrans,'FaceAlpha',0.2,'EdgeColor','black','EdgeAlpha',0.1);
%         clear ConHull;
%         clear ConVol;
%         clear denTree;
%         hold on;
%         treeCol(i,:) = ctrans;
%         scatter3(CellSoma(i,1),CellSoma(i,2),CellSoma(i,3),200, 'Marker','o', 'MarkerFaceColor',ctrans, 'MarkerEdgeColor','k','LineWidth', 2);
%         daspect([1 1 1]);
%         axis vis3d;
%         axis([20000 140000  60000 250000 -60000 0]);
%         plot( [20000, 40000], [70000, 70000],'-k' ) % insert 20um sclaebar
%         box off;
%         XColor = [1,1,1]; YColor = [1,1,1];
%         set (gca, 'XTick',[], 'YTick',[],'ZTick', [], 'Ydir','reverse');
%         view([-180,90]);
%     elseif ismember(cellIDs{i}, cellIDsDbx) ==1
%         figure(3);
%         denTree = DendriticTree(allTrees{i},i,cellIDs,cdbx,false);
%         [ConHull ConVol] = convhull(denTree(:,1),denTree(:,2),denTree(:,3),'simplify', true);
%         DbxHull{i} = ConHull;
%         DbxVol = [DbxVol, ConVol];
%         htri1 = trimesh(ConHull,denTree(:,1),denTree(:,2),denTree(:,3),'faceColor',cdbx,'FaceAlpha',0.2,'EdgeColor','black','EdgeAlpha',0.1);
%         hold on;
%         treeCol(i,:) = cdbx;
%         scatter3(CellSoma(i,1),CellSoma(i,2),CellSoma(i,3),200, 'Marker','o', 'MarkerFaceColor',cdbx, 'MarkerEdgeColor','k','LineWidth', 2);
%         daspect([1 1 1]);
%         axis vis3d;
%         axis([20000 140000  60000 250000 -60000 0]);
%         plot( [20000, 40000], [70000, 70000],'-k' ) % insert 20um sclaebar
%         box off;
%         XColor = [1,1,1]; YColor = [1,1,1];
%         set (gca, 'XTick',[], 'YTick',[],'ZTick', [], 'Ydir','reverse');
%         view([-180,90]);
%         clear ConHull;
%         clear ConVol;
%         clear denTree;
%     else
%         figure(4);;
%         denTree = DendriticTree(allTrees{i},i,cellIDs,cbarhl,false);
%         [ConHull ConVol] = convhull(denTree(:,1),denTree(:,2),denTree(:,3),'simplify', true);
%         BarhlHull{i} = ConHull;
%         BarhlVol = [BarhlVol, ConVol];
%         htri1 = trimesh(ConHull,denTree(:,1),denTree(:,2),denTree(:,3),'faceColor',cbarhl,'FaceAlpha',0.2,'EdgeColor','black','EdgeAlpha',0.1);
%         hold on;
%         treeCol(i,:) = cbarhl;
%         scatter3(CellSoma(i,1),CellSoma(i,2),CellSoma(i,3),200, 'Marker','o', 'MarkerFaceColor',cbarhl, 'MarkerEdgeColor','k','LineWidth', 2);
%         daspect([1 1 1]);
%         axis vis3d;
%         axis([20000 140000  60000 250000 -60000 0]);
%         plot( [20000, 40000], [70000, 70000],'-k' ) % insert 20um sclaebar
%         box off;
%         XColor = [1,1,1]; YColor = [1,1,1];
%         set (gca, 'XTick',[], 'YTick',[],'ZTick', [], 'Ydir','reverse');
%         view([-180,90]);
%         clear ConHull;
%         clear ConVol;
%         clear denTree;
%     end
%     pause(3);
% end

%%

figure();

AllConVol = [AlxVol, TransVol, DbxVol, BarhlVol];
AllColor = [repmat(calx,size(AlxVol,2),1);repmat(ctrans,size(TransVol,2),1); repmat(cdbx,size(DbxVol,2),1); repmat(cbarhl,size(BarhlVol,2),1)];

plot(repmat(1,1,length(AlxVol)), AlxVol./1e9, 'o','MarkerFaceColor',calx,'MarkerEdgeColor','k','MarkerSize',25);
hold on;
plot(repmat(2,1,length(TransVol)), TransVol./1e9, 'o', 'MarkerFaceColor',ctrans,'MarkerEdgeColor','k','MarkerSize',25);
plot(repmat(3,1,length(DbxVol)), DbxVol./1e9, 'o', 'MarkerFaceColor',cdbx,'MarkerEdgeColor','k','MarkerSize',25);
plot(repmat(4,1,length(BarhlVol)), BarhlVol./1e9, 'o', 'MarkerFaceColor',cbarhl,'MarkerEdgeColor','k','MarkerSize',25);


set(gca,'XTick', [1:4],'XTickLabel', {'group1'; 'group2'; 'group3'; 'group4'}, 'FontName', 'Arial', 'FontSize', 20);
set(gca, 'XLim', [0.5 4.5], 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
set(gcf,'color','w');
ylabel('Convexhull volume (\mum^3)', 'FontName', 'Arial', 'FontSize', 40);
axis square;
box off;
hold off;

%sum(AlxVol)/ (112*57*220*1e9) % total vol calculations

%%
figure();

histogram(AllConVol./1e9, 'BinWidth', min(AllConVol)/1e9);
ylabel('Count', 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
xlabel('Dendritic convex hull volume \mum^3', 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
set(gca, 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
axis square;
box off;

figure();
[y,I] = sort(AllConVol);
for i = 1:numel(cellIDs)
    plot(i, y(i)/1e9,'o','MarkerFaceColor',AllColor(I(i),:), 'MarkerEdgeColor','k', 'MarkerSize',35);
    hold on;
end
xlabel('cell #', 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
ylabel('Dendritic convex hull volume \mum^3', 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
set(gca, 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
axis square;



