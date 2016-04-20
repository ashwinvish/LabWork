peakAlx = [];
for i = 1:numel(cellIDs) 
    if ismember (cellIDs{i}, cellIDsAlx)==1
        [peak] = DendriticTree(allTrees{i},i,cellIDs,calx,true);
        peakAlx=[peakAlx,abs(CellSoma(i,3)/1000-peak)];
    else
        continue;
    end
end
clear peak;

figure;
peakTrans = [];
for i = 1:numel(cellIDs)
    if ismember (cellIDs{i}, cellIDsTrans)==1
       [peak] = DendriticTree(allTrees{i},i,cellIDs,ctrans,true);
       peakTrans = [peakTrans,abs(CellSoma(i,3)/1000 - peak)];
    else
        continue;
    end
end
clear peak;


figure();
peakDbx = [];
for i = 1:numel(cellIDs)
    if ismember (cellIDs{i}, cellIDsDbx)==1
        peak = DendriticTree(allTrees{i},i,cellIDs,cdbx,true);
        peakDbx = [peakDbx,abs(CellSoma(i,3)/1000-peak)];
    else
        continue;
    end
end
clear peak;

figure();
peakBarhl = [];
for i = 1:numel(cellIDs)
    if ismember (cellIDs{i}, cellIDsL)==1
        peak = DendriticTree(allTrees{i},i,cellIDs,cbarhl,true);
        peakBarhl = [peakBarhl, abs(CellSoma(i,3)/1000-peak)];
    else
        continue;
    end
end
clear peak;

figure;

% allPeaks = [peakAlx';peakTrans';peakDbx';peakBarhl'];
% groups = [ones(length(peakAlx),1); 2*ones(length(peakTrans),1); 3*ones(length(peakDbx),1); 4*ones(length(peakBarhl),1)];
% boxplot(allPeaks, groups);
% hold on;
% plot(ones(1,size(peakAlx,2)), peakAlx, 'o', 'MarkerFaceColor',calx, 'MarkerEdgeColor', 'k', 'MarkerSize', 35);
% plot(2*ones(1,size(peakTrans,2)), peakTrans, 'o','MarkerFaceColor',ctrans, 'MarkerEdgeColor', 'k', 'MarkerSize', 35);
% plot(3*ones(1,size(peakDbx,2)), peakDbx, 'o', 'MarkerFaceColor',cdbx, 'MarkerEdgeColor', 'k', 'MarkerSize', 35);
% plot(4*ones(1,size(peakBarhl,2)), peakBarhl, 'o', 'MarkerFaceColor',cbarhl, 'MarkerEdgeColor', 'k', 'MarkerSize', 35);
% 
% set(gca,'XLim',[0,5],'XTick',[1,2,3,4],'XTickLabel',{'group1', 'group2', 'group3', 'group4'}, 'YLim',[0 60],'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
% ylabel('Peak arbor depth from soma in \mum', 'FontName', 'Arial', 'FontSize', 40);
% set(gcf,'color','w');
% axis square;
% box off;





