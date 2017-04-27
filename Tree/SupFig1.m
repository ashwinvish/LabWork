% plot time constatants by spatial location

%Mediolateral
figure();
subplot(1,3,1);
scatter( CellSoma(:,1)/1000, log2(tau),750, CellColor, 'filled','MarkerEdgeColor','k');
xlabel('Medial --> Lateral', 'FontName', 'Arial', 'FontSize', 40);
ylabel('log(\tau)','FontName', 'Arial', 'FontSize', 40);
set(gca,'YLim', [0, max(log2(tau))], 'Xlim', [min(CellSoma(:,1)/1000), max(CellSoma(:,1)/1000)], 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
set(gcf,'color','w');
box off;
axis square;
[R1,P1] = corrcoef( CellSoma(:,1), log2(tau));

%RostroCaudal
subplot(1,3,2);
scatter( CellSoma(:,2)/1000, log2(tau),750, CellColor, 'filled','MarkerEdgeColor','k');
xlabel('Caudal --> Rostral', 'FontName', 'Arial', 'FontSize', 40);
ylabel('log(\tau)','FontName', 'Arial', 'FontSize', 40);
set(gca,'YLim', [0, max(log2(tau))], 'Xlim', [min(CellSoma(:,2)/1000), max(CellSoma(:,2)/1000)], 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
set(gcf,'color','w');
box off;
axis square;
[R2,P2] = corrcoef( CellSoma(:,2), log2(tau));

% Dorsoventral
subplot(1,3,3);
scatter( CellSoma(:,3)/1000, log2(tau),750, CellColor, 'filled','MarkerEdgeColor','k');
xlabel('Dorsal --> Ventral', 'FontName', 'Arial', 'FontSize', 40);
ylabel('log(\tau)','FontName', 'Arial', 'FontSize', 40);
set(gca,'YLim', [0, max(log2(tau))], 'Xlim', [min(CellSoma(:,3)/1000), max(CellSoma(:,3)/1000)], 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
set(gcf,'color','w');
box off;
axis square;
[R3,P3] = corrcoef( CellSoma(:,3), log2(tau));

% sorted time constants







