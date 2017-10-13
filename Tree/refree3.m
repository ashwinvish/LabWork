% plot time constants and synapse numbers

AllPostLength =  cellfun(@length, allPost);
CellGroupColor;
% plot number of postsynapse vs time constants

scatter(AllPostLength,rho,500,CellColor,'filled','MarkerEdgeColor','k');
set(gca,'XLim',[0,200],'YLim',[0, 7],'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
ylabel('log(\tau)', 'FontName', 'Arial', 'FontSize', 40);
xlabel('Number of postsynaptic sites','FontName', 'Arial', 'FontSize', 40);
%set(gcf,'color','w');
axis square;
%box off;

[R1,P1] = corrcoef(AllPostLength,log2(tau));

SynDensity = AllPostLength./denLength;
figure;
scatter(SynDensity*1000,rho,500,CellColor,'filled','MarkerEdgeColor','k');
set(gca,'XLim',[0,1],'YLim',[0, 7],'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
ylabel('log(\tau)', 'FontName', 'Arial', 'FontSize', 40);
xlabel('Synaptic Density','FontName', 'Arial', 'FontSize', 40);
%set(gcf,'color','w');
axis square;
%box off;

[R2,P2] = corrcoef(SynDensity,log2(tau));

figure;
scatter3(SynDensity*1000,AllPostLength,rho,500,CellColor,'filled','MarkerEdgeColor','k');
% set(gca,'XLim',[0,1],'YLim',[0, 7],'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
% ylabel('log(\tau)', 'FontName', 'Arial', 'FontSize', 40);


