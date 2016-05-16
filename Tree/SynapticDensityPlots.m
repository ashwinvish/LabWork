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




