AlxSomaNeurites =[];
TransSomaNeurites = [];
DbxSomaNeurites = [];
BarhlSomaNeurites= [];

for i = 1:numel(cellIDs)
    if ismember(cellIDs{i},cellIDsAlx) ==1
        AlxSomaNeurites = [AlxSomaNeurites,length(allTrees{i}{1,1}{2})];
        plot(i,length(allTrees{i}{1,1}{2}), 'o','MarkerFaceColor', calx, 'MarkerEdgeColor', 'k', 'MarkerSize', 25);
    elseif ismember(cellIDs{i},cellIDsTrans) ==1
        TransSomaNeurites = [TransSomaNeurites,length(allTrees{i}{1,1}{2})];
        plot(i,length(allTrees{i}{1,1}{2}), 'o','MarkerFaceColor', ctrans, 'MarkerEdgeColor', 'k', 'MarkerSize', 25);
    elseif ismember(cellIDs{i},cellIDsDbx) ==1
        DbxSomaNeurites = [DbxSomaNeurites,length(allTrees{i}{1,1}{2})];
        plot(i,length(allTrees{i}{1,1}{2}), 'o','MarkerFaceColor', cdbx, 'MarkerEdgeColor', 'k', 'MarkerSize', 25);
    else
        BarhlSomaNeurites = [BarhlSomaNeurites,length(allTrees{i}{1,1}{2})];
        plot(i,length(allTrees{i}{1,1}{2}), 'o','MarkerFaceColor', cbarhl, 'MarkerEdgeColor', 'k', 'MarkerSize', 25);
    end
    hold on;
end

xlabel('Cell #','FontName', 'Arial', 'FontSize', 40);
ylabel('Neurite from somata','FontName', 'Arial', 'FontSize', 40);
set(gca, 'XLim',[0 22], 'YLim',[0 10],'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
axis square; box off;

SomaNeurites = [AlxSomaNeurites,TransSomaNeurites,DbxSomaNeurites,BarhlSomaNeurites ]
SomaMeanNeurites = [mean(AlxSomaNeurites);mean(TransSomaNeurites);mean(DbxSomaNeurites);mean(BarhlSomaNeurites)];
SomaNeuriteError = [std(AlxSomaNeurites);std(TransSomaNeurites);std(DbxSomaNeurites);std(BarhlSomaNeurites)];


