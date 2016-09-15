%% New Plots

% RC vs Tau
figure();
plot(CellSoma(:,2)/1000, log2(tau), 'ko','markersize',25,'LineWidth', 2);
xlabel('Rostral --> caudal (\mum)');
ylabel('log(tau)');
set(gca,'XLim',[min(CellSoma(:,2)/1000) max(CellSoma(:,2)/1000)], 'FontName', 'Arial','FontSize',40, 'LineWidth', 2)
axis square;
box off;


% ML vs Tau

figure();
plot(CellSoma(:,1)/1000, log2(tau), 'ko','markersize',25,'LineWidth', 2);
xlabel('Medial --> Lateral (\mum)');
ylabel('log(tau)');
set(gca,'XLim',[min(CellSoma(:,1)/1000), max(CellSoma(:,1)/1000)],'FontName', 'Arial','FontSize',40, 'LineWidth', 2)
axis square;
box off;