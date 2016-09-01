for i = 1:numel(cellIDs)
    figure(3);
    plot(i,rho(i),'o','MarkerFaceColor', MarkerColorMap(i,:), 'MarkerEdgeColor','k','MarkerSize',35);
    hold on;
end

xlabel('cell#','FontName', 'Arial', 'FontSize', 40);
ylabel('Normalized perisistence','FontName', 'Arial', 'FontSize', 40);
set(gca, 'XLim',[0 22], 'YLim',[0 1],'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
axis square; box off;