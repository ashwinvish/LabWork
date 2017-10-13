%%
SomaDiam = [2.90775
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

run CellGroupColor.m
MarkerColorMap = CellColor;
[y,I] = sort(rho);

figure()
subplot(1,3,1);
scatter(1:22,y,750,CellColor(I,:),'filled','MarkerEdgeColor','k');
xlabel('Soma diameter (\mum)','FontName', 'Arial', 'FontSize', 40);
ylabel('Normalized perisistence \rho','FontName', 'Arial', 'FontSize', 40);
set(gca,'XLIm', [0,22], 'YLim',[0 7],'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
axis square; box off;

subplot(1,3,2);
scatter(CellSoma(:,3)/1000,SomaDiam,750,CellColor,'filled','MarkerEdgeColor','k');
xlabel('Soma depth (\mum)','FontName', 'Arial', 'FontSize', 40);
ylabel('Soma diameter (\mum)','FontName', 'Arial', 'FontSize', 40);
set(gca, 'YLim',[0 max(SomaDiam)],'XLim',[0 max(CellSoma(:,3))/1000],'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
axis square; box off;

subplot(1,3,3);
scatter(SomaDiam,rho, 750, MarkerColorMap,'filled', 'MarkerEdgeColor','k');
xlabel('SomaDiameter','FontName', 'Arial', 'FontSize', 40);
ylabel('Log time constant','FontName', 'Arial', 'FontSize', 40);
set(gca, 'XLim',[min(SomaDiam) max(SomaDiam)], 'YLim',[0 7],'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
axis square; box off

hold on;

X = [ones(length(SomaDiam),1), SomaDiam];
b = X\rho';
yCalc2 = X*b;
plot( SomaDiam,yCalc2,'-r');

R2 = 1 - sum((rho' - yCalc2).^2)/sum((rho'-mean(rho')).^2);
text(1,1, R2);

%%

AlxSomaDia = [];
TransSomaDia = [];
DbxSomaDia = [];
BarhlSomaDia = [];

for i = 1:numel(cellIDs)
    if ismember(cellIDs{i}, cellIDsAlx) ==1
        AlxSomaDia = [AlxSomaDia; SomaDiam(i)];
        BarCMap(i,:) = calx;
    elseif ismember( cellIDs{i}, cellIDsTrans) ==1
        TransSomaDia = [TransSomaDia; SomaDiam(i)];
        BarCMap(i,:) = ctrans;
    elseif ismember( cellIDs{i}, cellIDsDbx) ==1
        DbxSomaDia = [DbxSomaDia; SomaDiam(i)];
        BarCMap(i,:) = cdbx;
    else
        BarhlSomaDia = [BarhlSomaDia;SomaDiam(i)];
        BarCMap(i,:) = cbarhl;
    end
end