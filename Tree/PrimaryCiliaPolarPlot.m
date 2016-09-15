PrimaryCilia = [allTrees{1}{1,2}{1,3};
    allTrees{2}{1,2}{1,3};
    allTrees{3}{1,73}{1,3};
    allTrees{4}{1,2}{1,3};
    allTrees{5}{1,2}{1,3};
    allTrees{6}{1,3}{1,3};
    allTrees{7}{1,2}{1,3};
    allTrees{8}{1,8}{1,3};
    allTrees{9}{1,2}{1,3};
    allTrees{10}{1,2}{1,3};
    allTrees{11}{1,3}{1,3};
    allTrees{12}{1,3}{1,3};
    allTrees{13}{1,2}{1,3};
    allTrees{13}{1,3}{1,3};
    allTrees{15}{1,3}{1,3};
    allTrees{16}{1,2}{1,3};
    [5*19711, 5* 23075, 45*436];
    allTrees{18}{1,5}{1,3};
    allTrees{19}{1,3}{1,3};
    allTrees{20}{1,4}{1,3};
    allTrees{21}{1,3}{1,3};
    allTrees{22}{1,3}{1,3}];

for i = 1:numel(cellIDs)
    if ismember(cellIDs{i}, cellIDsAlx) ==1
        subplot(221)
        h1 = compass(PrimaryCilia(i,1)/1000-CellSoma(i,1)/1000, PrimaryCilia(i,2)/1000-CellSoma(i,2)/1000);
        h1.LineWidth = 2;
        h1.Color = calx;
        hold on;
    end
end
% axis([min(h1.XData) max(h1.XData), min(h1.YData), max(h1.YData)])

for i = 1:numel(cellIDs)
    if ismember(cellIDs{i}, cellIDsTrans) ==1
        subplot(222)
        h2 = compass(PrimaryCilia(i,1)/1000-CellSoma(i,1)/1000, PrimaryCilia(i,2)/1000-CellSoma(i,2)/1000);
        h2.LineWidth = 2;
        h2.Color = ctrans;
        hold on;
    end
end
% axis([min(h2.XData) max(h2.XData), min(h2.YData), max(h2.YData)])


for i = 1:numel(cellIDs)

    if ismember(cellIDs{i}, cellIDsDbx) ==1
        subplot(223)
        h3 = compass(PrimaryCilia(i,1)/1000-CellSoma(i,1)/1000, PrimaryCilia(i,2)/1000-CellSoma(i,2)/1000);
        h3.LineWidth = 2;
        h3.Color = cdbx;
        hold on;
    end
end
axis([min(h3.XData) max(h3.XData), min(h3.YData), max(h3.YData)])


for i = 1:numel(cellIDs)
    if ismember(cellIDs{i}, cellIDsL) ==1
        subplot(224)
        h4 = compass(PrimaryCilia(i,1)/1000-CellSoma(i,1)/1000, PrimaryCilia(i,2)/1000-CellSoma(i,2)/1000);
        h4.LineWidth = 2;
        h4.Color = cbarhl;
        hold on;
    end
 
end
%     axis([min(h4.XData) max(h4.XData), min(h4.YData), max(h4.YData)])





