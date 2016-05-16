function [  ] = GroupTaper( groupAxonPathlength, groupDenPathlength, groupAxnDia, groupDenDia )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

groupAxonPathlength(groupAxonPathlength == 0) = NaN;
groupDenPathlength(groupDenPathlength == 0) = NaN;
groupDenDia(groupDenDia == 0) = NaN;
groupAxnDia(groupAxnDia == 0) = NaN;

steps = [0:0.1:1];

for i = 1:1:10
    id1 = find(groupAxonPathlength > steps(i) & groupAxonPathlength < steps(i+1));
    id2 = find(groupDenPathlength> steps(i)  &  groupDenPathlength< steps(i+1));
    
    groupAxnPathLengthMean(i) = nanmean(groupAxonPathlength(id1));
    groupAxnPathLengthSD(i) = nanstd(groupAxonPathlength(id1));
    
    groupAxonDiaMean(i) = nanmean(groupAxnDia(id1));
    groupAxonDiaSD(i) = nanstd(groupAxnDia(id1));
    
    groupDenPlengthMean(i) = nanmean(groupDenPathlength(id2));
    groupDenPlengthSD(i) = nanstd(groupDenPathlength(id2));
    
    groupDenDiaMean(i) = nanmean(groupDenDia(id2));
    groupDenDiaSD(i) = nanstd(groupDenDia(id2));
    
end

errorbar(groupAxnPathLengthMean,groupAxonDiaMean, groupAxonDiaSD, '-o','Color',[0,0.8,0], 'LineWidth',2);
hold on;
errorbar(groupDenPlengthMean,groupDenDiaMean, groupDenDiaSD,'-o', 'Color',[.9,0,0],'LineWidth',2);
xlabel('Norm. pathlength', 'FontName', 'Arial', 'FontSize', 40);
ylabel('Dia. (nm)','FontName', 'Arial', 'FontSize', 40);
%legend({'Axon','Dendrite'},'FontName', 'Arial', 'FontSize', 40);
set(gca,'XLim',[0 1], 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
set(gcf,'color','w');
box off;
hold off;
axis square;
end

