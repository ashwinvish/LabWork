
%%
AlxSomaDia = [];
TransSomaDia = [];
DbxSomaDia = [];
BarhlSomaDia = [];

for i = 1:numel(cellIDs)
    if ismember(cellIDs{i}, cellIDsAlx) ==1
        AlxSomaDia = [AlxSomaDia; SomaDiameters(i)];
        BarCMap(i,:) = calx;
    elseif ismember( cellIDs{i}, cellIDsTrans) ==1
        TransSomaDia = [TransSomaDia; SomaDiameters(i)];
        BarCMap(i,:) = ctrans;
    elseif ismember( cellIDs{i}, cellIDsDbx) ==1
        DbxSomaDia = [DbxSomaDia; SomaDiameters(i)];
        BarCMap(i,:) = cdbx;
    else 
        BarhlSomaDia = [BarhlSomaDia;SomaDiameters(i)];
        BarCMap(i,:) = cbarhl;
    end
end

%%
SomaDiameters = [2.90775
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
[y,I] = sort(rho);

figure()


for i = 1:numel(cellIDs)
    subplot(1,3,1)
    hold on;
    plot(i, rho(I(i)), 'o','MarkerFaceColor', MarkerColorMap(I(i),:), 'MarkerEdgeColor','k','MarkerSize',35);
    
end

% plot(SomaDiameters, rho,'o','MarkerFaceColor', 'k', 'MarkerEdgeColor','w','MarkerSize',35);
% f = ezfit('exp');
% showfit(f, 'fitcolor', 'red', 'fitlinewidth',2);




% P = polyfit(SomaDiameter,rho',1)
% x1 = min(SomaDiameter):0.5:7;
% y1= polyval(P,x1);
% plot(x1,y1, '--k');
% hold off;

xlabel('Soma diameter (\mum)','FontName', 'Arial', 'FontSize', 40);
ylabel('Normalized perisistence \rho','FontName', 'Arial', 'FontSize', 40);
set(gca,'XLIm', [0,22], 'YLim',[0 7],'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
axis square; box off;


for i = 1:numel(cellIDs)
    
subplot(1,3,2);
    plot(CellSoma(i,3)/1000,SomaDiameters(i),'o','MarkerFaceColor', MarkerColorMap(i,:), 'MarkerEdgeColor','k','MarkerSize',35);
    hold on;
end

xlabel('Soma depth (\mum)','FontName', 'Arial', 'FontSize', 40);
ylabel('Soma diameter (\mum)','FontName', 'Arial', 'FontSize', 40);
set(gca, 'YLim',[0 max(SomaDiameters)],'XLim',[0 max(CellSoma(:,3))/1000],'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
axis square; box off;

subplot(1,3,3);
scatter(SomaDiameters,rho, 1000, MarkerColorMap,'filled', 'MarkerEdgeColor','k');
xlabel('SomaDiameter','FontName', 'Arial', 'FontSize', 40);
ylabel('Log time constant','FontName', 'Arial', 'FontSize', 40);
set(gca, 'XLim',[min(SomaDiameters) max(SomaDiameters)], 'YLim',[0 7],'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
axis square; box off

hold on;

X = [ones(length(SomaDiameters),1), SomaDiameters];
b = X\rho';
yCalc2 = X*b;
plot( SomaDiameters,yCalc2,'-r');

R2 = 1 - sum((rho' - yCalc2).^2)/sum((rho'-mean(rho')).^2);
text(1,1, R2);


        



