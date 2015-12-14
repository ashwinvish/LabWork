hold on;
for i = 1:length(cellIDs)
    if ismember(cellIDs(i),cellIDsAlx) == 1
        BarCMap = calx;
        scatter3(CellSoma(i,1),CellSoma(i,2),CellSoma(i,3),500,'MarkerFaceColor',BarCMap, 'MarkerEdgeColor', 'k');
    elseif ismember(cellIDs(i),cellIDsDbx) == 1
        BarCMap = cdbx;
        scatter3(CellSoma(i,1),CellSoma(i,2),CellSoma(i,3),500,'MarkerFaceColor',BarCMap, 'MarkerEdgeColor', 'k');
    else
        BarCMap= cbarhl;
        scatter3(CellSoma(i,1),CellSoma(i,2),CellSoma(i,3),500,'MarkerFaceColor',BarCMap, 'MarkerEdgeColor', 'k');
    end
end
box on;
%axis([ 20000 140000 60000 250000 -60000 0]);
%plot( [20000, 40000], [70000, 70000],'-k' ) % insert 20um sclaebar
%scatter3(MauthnerCell(1,1),MauthnerCell(1,2),MauthnerCell(1,3), 500,'MarkerFaceColor','k', 'MarkerEdgeColor', 'k'); % location of the Mauther Cell center    
daspect([1 1 1]); % make aspect ratio [1 1 1]
%set (gca,'Ydir','reverse');
set (gca,'XTick',[], 'YTick',[],'ZTick', [], 'Ydir','reverse');
%view([-180,90]); % xy view