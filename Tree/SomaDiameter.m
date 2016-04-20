SomaDiameter = [2.90775
4.8832
5.66775
4.18275
4.4138
3.5702
4.15525
4.40295
5.23905
5.0998
6.00455
4.65965
4.738
4.178
4.97925
4.0166
4.881
4.85725
4.699
4.18465
3.83845
4.68705];

for i = 1:numel(cellIDs)
    plot(SomaDiameter(i), rho(i), 'o','MarkerFaceColor', treeCol(i,:), 'MarkerEdgeColor','k','MarkerSize',25);
    hold on;
end
xlabel('Soma diameter (\mum)','FontName', 'Arial', 'FontSize', 40);
ylabel('Normalized perisistence \rho','FontName', 'Arial', 'FontSize', 40);
set(gca, 'XLim',[0 max(SomaDiameter)], 'YLim',[0 1],'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
axis square; box off;