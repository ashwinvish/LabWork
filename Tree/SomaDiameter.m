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

MarkerColorMap = BarCMap;

for i = 1:numel(cellIDs)
    figure(1);
    plot(SomaDiameter(i), rho(i), 'o','MarkerFaceColor', MarkerColorMap(i,:), 'MarkerEdgeColor','k','MarkerSize',35);
    hold on;
end

xlabel('Soma diameter (\mum)','FontName', 'Arial', 'FontSize', 40);
ylabel('Normalized perisistence \rho','FontName', 'Arial', 'FontSize', 40);
set(gca, 'XLim',[0 max(SomaDiameter)], 'YLim',[0 1],'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
axis square; box off;

for i = 1:numel(cellIDs)
    figure(2)
    plot(CellSoma(i,3)/1000,SomaDiameter(i),'o','MarkerFaceColor', MarkerColorMap(i,:), 'MarkerEdgeColor','k','MarkerSize',35);
    hold on;
end

xlabel('Soma depth (\mum)','FontName', 'Arial', 'FontSize', 40);
ylabel('Soma diameter (\mum)','FontName', 'Arial', 'FontSize', 40);
set(gca, 'YLim',[0 max(SomaDiameter)],'XLim',[0 max(CellSoma(:,3))/1000],'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
axis square; box off;


% X = [ones(length(SomaDiameter),1), SomaDiameter];
% b = X\rho';
% yCalc2 = X*b;
% plot( SomaDiameter,yCalc2,'-r');
% 
% P = polyfit(SomaDiameter,rho',1)
% x1 = min(SomaDiameter):0.5:7;
% y1= polyval(P,x1);
% plot(x1,y1, '--k');
% hold off;



hold off;




for i = 1:numel(cellIDs)
    figure(3);
    plot(i,rho(i),'o','MarkerFaceColor', MarkerColorMap(i,:), 'MarkerEdgeColor','k','MarkerSize',35);
    hold on;
    
end
xlabel('cell#','FontName', 'Arial', 'FontSize', 40);
ylabel('Normalized perisistence','FontName', 'Arial', 'FontSize', 40);
set(gca, 'XLim',[0 22], 'YLim',[0 1],'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
axis square; box off;

figure(4);

scatter(1:numel(cellIDs), sort(rho),1000, MarkerColorMap,'filled', 'MarkerEdgeColor','k');
xlabel('cell#','FontName', 'Arial', 'FontSize', 40);
ylabel('Perisistence','FontName', 'Arial', 'FontSize', 40);
set(gca, 'XLim',[0 22], 'YLim',[0 1],'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
axis square; box off





%%
AlxSomaDia = [];
TransSomaDia = [];
DbxSomaDia = [];
BarhlSomaDia = [];

for i = 1:numel(cellIDs)
    if ismember(cellIDs{i}, cellIDsAlx) ==1
        AlxSomaDia = [AlxSomaDia; SomaDiameter(i)];
        BarCMap(i,:) = calx;
    elseif ismember( cellIDs{i}, cellIDsTrans) ==1
        TransSomaDia = [TransSomaDia; SomaDiameter(i)];
        BarCMap(i,:) = ctrans;
    elseif ismember( cellIDs{i}, cellIDsDbx) ==1
        DbxSomaDia = [DbxSomaDia; SomaDiameter(i)];
        BarCMap(i,:) = cdbx;
    else 
        BarhlSomaDia = [BarhlSomaDia;SomaDiameter(i)];
        BarCMap(i,:) = cbarhl;
    end
end

        



