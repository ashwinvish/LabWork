AlxSomaTime = [];
TransSomaTime = [];
DbxSomaTime  = [];
BarhlSomaTime = [];
CellColor = [];

for i = 1:numel(cellIDs)
    if ismember(cellIDs{i},cellIDsAlx) ==1
        AlxSomaTime = [AlxSomaTime; tau(i)];
    elseif ismember(cellIDs{i}, cellIDsTrans) ==1
        TransSomaTime = [TransSomaTime; tau(i)];
    elseif ismember(cellIDs{i}, cellIDsDbx) ==1
        DbxSomaTime = [DbxSomaTime; tau(i)];
    else
        BarhlSomaTime = [BarhlSomaTime; tau(i)];
    end
end


h = boxplot([log2(AlxSomaTime);log2(TransSomaTime);log2(DbxSomaTime);log2(BarhlSomaTime)],[ones(size(AlxSomaTime,1),1); 2*ones(size(TransSomaTime,1),1); ...
     3*ones(size(DbxSomaTime,1),1); 4*ones(size(BarhlSomaTime,1),1)],...
'Notch','off', 'Symbol', 'ko','Colors',[calx;ctrans;cdbx;cbarhl],'OutlierSize', 25);
set(h,{'linew'},{4});
h1 = findobj(gca,'tag','Median');
set(h1,'Color','k');
ylabel('log(time)', 'FontName', 'Arial', 'FontSize', 40);
set(gca,'YLim',[0 7],'FontName', 'Arial', 'FontSize', 40, 'LineWidth',4);
box off;

hold on;
plot([ones(size(AlxSomaTime,1),1)', 2*ones(size(TransSomaTime,1),1)',3*ones(size(DbxSomaTime,1),1)', 4*ones(size(BarhlSomaTime,1),1)' ]...
    ,[log2(AlxSomaTime);log2(TransSomaTime);log2(DbxSomaTime);log2(BarhlSomaTime)]',...
    'o','MarkerSize', 25, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'k', 'LineWidth', 4);
set(findobj(gcf,'LineStyle','--'),'LineStyle','-');


hold off;
axis square;