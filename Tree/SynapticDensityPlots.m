for i = 1:length(cellIDs)
    if ismember(cellIDs{i},cellIDsAlx) == 1
        BarCMap(i,:) = calx;
    elseif ismember(cellIDs{i},cellIDsTrans) == 1
        BarCMap(i,:) = ctrans;
    elseif ismember(cellIDs{i},cellIDsDbx) == 1
        BarCMap(i,:) = cdbx;
    else
        BarCMap(i,:)= cbarhl;
    end
end



[m,n] = sort(denLength);
[y,I] = sort(cellfun(@length,allPost));
subplot(1,2,1);
scatter(m/1000,y(n),750,BarCMap(n,:),'filled','MarkerEdgeColor','k');
xlabel('Dendrite length', 'FontName', 'Arial', 'FontSize', 40);
ylabel('# postsynaptic site','FontName', 'Arial', 'FontSize', 40);
set(gca,'XLim', [0, max(m/1000)], 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
set(gcf,'color','w');
box off;
axis square;

subplot(1,2,2);
scatter(1:22,y(n)./(m/1000),750,BarCMap(n,:),'filled','MarkerEdgeColor','k');
xlabel('cell #', 'FontName', 'Arial', 'FontSize', 40);
ylabel('density (#sites/den. length)','FontName', 'Arial', 'FontSize', 40);
set(gca,'XLim', [0, 22], 'YLim', [0 1], 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
set(gcf,'color','w');
box off;
axis square;



figure();

temp = cell2mat(allRawLength);
subplot(1,2,1);
scatter(temp(n)/1000,y(n),750,BarCMap(n,:),'filled','MarkerEdgeColor','k');
xlabel('Raw length', 'FontName', 'Arial', 'FontSize', 40);
ylabel('# postsynaptic site','FontName', 'Arial', 'FontSize', 40);
set(gca,'XLim',[0 max(temp(n)/1000)], 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
set(gcf,'color','w');
box off;
axis square;

subplot(1,2,2);
scatter(1:22,y(n)./(temp(n)/1000),750,BarCMap(n,:),'filled','MarkerEdgeColor','k');
xlabel('cell #', 'FontName', 'Arial', 'FontSize', 40);
ylabel('density (#sites/raw length)','FontName', 'Arial', 'FontSize', 40);
set(gca,'XLim', [0, 22], 'YLim', [0 1], 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
set(gcf,'color','w');
box off;
axis square;


%%
AlxPostSynapseDensity = [];
TransPostSynapseDensity = [];
DbxPostSynapseDensity = [];
BarhlPostSynapseDensity = [];


for i = 1:numel(cellIDs)
    if ismember(cellIDs{i}, cellIDsAlx) ==1
        AlxPostSynapseDensity = [AlxPostSynapseDensity; denLength(i)/ length(allPost{i})];
    elseif ismember(cellIDs{i}, cellIDsTrans) ==1
        TransPostSynapseDensity = [TransPostSynapseDensity;  denLength(i)/length(allPost{i})];
    elseif ismember(cellIDs{i}, cellIDsDbx) ==1
        DbxPostSynapseDensity = [DbxPostSynapseDensity; denLength(i)/length(allPost{i})];
    else
        BarhlPostSynapseDensity = [ BarhlPostSynapseDensity; denLength(i)/length(allPost{i})];
    end
end

meanDensity = [mean(AlxPostSynapseDensity/1000); mean(TransPostSynapseDensity/1000); mean(DbxPostSynapseDensity/1000); mean(BarhlPostSynapseDensity/1000)];
StdDensity = [std(AlxPostSynapseDensity/1000); std(TransPostSynapseDensity/1000); std(DbxPostSynapseDensity/1000); std(BarhlPostSynapseDensity/1000)];


plot(ones(length(AlxPostSynapseDensity),1), AlxPostSynapseDensity/1000, 'o','MarkerSize', 35, 'MarkerFaceColor', calx, 'MarkerEdgeColor', 'k');
hold on;
plot(2*ones(length(TransPostSynapseDensity),1), TransPostSynapseDensity/1000, 'o','MarkerSize', 35, 'MarkerFaceColor', ctrans, 'MarkerEdgeColor', 'k');
plot(3*ones(length(DbxPostSynapseDensity),1), DbxPostSynapseDensity/1000, 'o','MarkerSize', 35, 'MarkerFaceColor', cdbx, 'MarkerEdgeColor', 'k');
plot(4*ones(length(BarhlPostSynapseDensity),1), BarhlPostSynapseDensity/1000, 'o','MarkerSize', 35, 'MarkerFaceColor', cbarhl, 'MarkerEdgeColor', 'k');
plot([1,2,3,4],meanDensity, 'o','MarkerSize', 35, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
plot([1:4;1:4], [meanDensity'-StdDensity';meanDensity'+StdDensity'], 'Color','k','LineWidth',2);


axis square;
box off;
set(gca, 'XTick', [1:4],'XTickLabel', {'Ipsi'; 'Ipsi-Contra'; 'Contra'; 'Unknown'},'XLim', [0.5 4.5], 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
ylabel('#Synapses/\mum', 'FontName', 'Arial', 'FontSize', 40);
set(gcf,'color','w');







