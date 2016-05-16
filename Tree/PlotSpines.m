AlxSpines = [];
TransSpines = [];
DbxSpines = [];
BarhlSpines = [];

for i = 1:numel(cellIDs)
    if ismember(cellIDs{i}, cellIDsAlx) ==1
        AlxSpines = [AlxSpines, size(allSpine{i},1)];
    elseif ismember(cellIDs{i}, cellIDsTrans) ==1
        TransSpines = [TransSpines, size(allSpine{i},1)];
    elseif ismember(cellIDs{i}, cellIDsDbx) ==1
        DbxSpines = [DbxSpines, size(allSpine{i},1)];
    else
        BarhlSpines = [BarhlSpines, size(allSpine{i},1)];
    end
end

plot(ones(1,length(AlxSpines)),AlxSpines, 'o','MarkerFaceColor', calx, 'MarkerSize', 25, 'MarkerEdgeColor', 'k');
hold on;
plot(2*ones(1, length(TransSpines)), TransSpines, 'o','MarkerFaceColor', ctrans, 'MarkerSize', 25, 'MarkerEdgeColor', 'k');
plot(3*ones(1, length(DbxSpines)), DbxSpines,'o','MarkerFaceColor', cdbx, 'MarkerSize', 25, 'MarkerEdgeColor', 'k');
plot(4*ones(1, length(BarhlSpines)), BarhlSpines, 'o','MarkerFaceColor', cbarhl, 'MarkerSize', 25, 'MarkerEdgeColor', 'k');

set(gca, 'XLim', [0.5, 4.5],'XTick', [1,2,3,4], 'XTickLabel', {'Ipsi', 'Ipsi-Contra', 'Contra', 'Unknown'}, 'FontName', 'Arial', 'FontSize', 40, 'LineWidth', 2);
%xlabel({'Ipsi', 'Ipsi-Contra', 'Contra', 'Unknown'},  'FontName', 'Arial', 'FontSize', 40);
ylabel('Number of spinules',  'FontName', 'Arial', 'FontSize', 40);
box off;
axis square;



