%% New Plots

% RC vs Tau
figure();
subplot(1,3,1)
for i = 1:22
    plot(CellSoma(i,2)/1000, log2(tau(i)),'o', 'MarkerFaceColor', BarCMap(i,:),'MarkerEdgeColor', 'k','markersize',25,'LineWidth', 2);
        hold on;
end
xlabel('Rostral --> caudal (\mum)');
ylabel('log(tau)');
set(gca,'XLim',[min(CellSoma(:,2)/1000) max(CellSoma(:,2)/1000)],'YLim', [0,7], 'FontName', 'Arial','FontSize',40, 'LineWidth', 2)
axis square;
box off;


% ML vs Tau
subplot(1,3,2)
for i = 1:22
plot(CellSoma(i,1)/1000, log2(tau(i)),'o', 'MarkerFaceColor', BarCMap(i,:),'MarkerEdgeColor', 'k','markersize',25,'LineWidth', 2);
hold on;
end
xlabel('Medial --> Lateral (\mum)');
ylabel('log(tau)');
set(gca,'XLim',[min(CellSoma(:,1)/1000), max(CellSoma(:,1)/1000)],'YLim', [0,7],'FontName', 'Arial','FontSize',40, 'LineWidth', 2)
axis square;
box off;

% DV vs Tau
subplot(1,3,3)
for i = 1:22
plot(CellSoma(i,3)/1000, log2(tau(i)),'o', 'MarkerFaceColor', BarCMap(i,:),'MarkerEdgeColor', 'k','markersize',25,'LineWidth', 2);
hold on;
end
xlabel('Dorsal --> Ventral (\mum)');
ylabel('log(tau)');
set(gca,'XLim',[min(CellSoma(:,3)/1000), max(CellSoma(:,3)/1000)],'YLim', [0,7],'FontName', 'Arial','FontSize',40, 'LineWidth', 2)
axis square;
box off;

%%
% ML orginization of dendrites

[AlxXY, AlxYZ, AlxXZ] = DendriticTreeOfGroups(cellIDsAlx, allTrees, cellIDs, calx, false);
[TransXY, TransYZ, TransXZ] = DendriticTreeOfGroups(cellIDsTrans, allTrees, cellIDs, ctrans, false);
[DbxXY, DbxYZ, DbxXZ] = DendriticTreeOfGroups(cellIDsDbx, allTrees, cellIDs, cdbx, false);
[BarhlXY, BarhlYZ, BarhlXZ] = DendriticTreeOfGroups(cellIDsL, allTrees, cellIDs, cbarhl, false);

% normalize 
NormAlxXY = bsxfun(@rdivide, AlxXY,max(AlxXY,[],2));
NormTransXY = bsxfun(@rdivide, TransXY,max(TransXY,[],2));
NormDbxXY = bsxfun(@rdivide, DbxXY,max(DbxXY,[],2));
NormBarhlXY = bsxfun(@rdivide, BarhlXY,max(BarhlXY,[],2));

X= 20:5:140;

figure();

shadedErrorBar(X,NormAlxXY,{@nanmean, @nanstd},{'-','color',calx,'markerfacecolor',calx,'LineWidth',4},1);
hold on;
shadedErrorBar(X,NormTransXY,{@nanmean, @nanstd},{'-','color',ctrans,'markerfacecolor',ctrans,'LineWidth',4},1);
shadedErrorBar(X,NormDbxXY,{@nanmean, @nanstd},{'-','color',cdbx,'markerfacecolor',cdbx,'LineWidth',4},1);
shadedErrorBar(X,NormBarhlXY,{@nanmean, @nanstd},{'-','color',cbarhl,'markerfacecolor',cbarhl,'LineWidth',4},1);
set(gca,'XLim',[20 140],'XDir','reverse','XTickLabel',[0,20,40,60,80,100,120],'YTick',[0.5,1],'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2 );
box off;
axis square;

figure();

errorbar(X,nanmean(NormAlxXY,1), nanstd(NormAlxXY,1), '-o','Color', calx, 'MarkerEdgeColor', calx,'MarkerFaceColor', 'none', 'LineWidth', 2);
hold on;
%errorbar(X,nanmean(NormTransXY,1), nanstd(NormTransXY,1), '-o','Color', ctrans, 'MarkerEdgeColor', ctrans,'MarkerFaceColor', 'none');
errorbar(X,nanmean(NormDbxXY,1), nanstd(NormDbxXY,1), '-o','Color', cdbx, 'MarkerEdgeColor', cdbx,'MarkerFaceColor', 'none', 'LineWidth', 2);
%errorbar(X,nanmean(NormBarhlXY,1), nanstd(NormBarhlXY,1), '-o','Color', cbarhl, 'MarkerEdgeColor', cbarhl,'MarkerFaceColor', 'none');
set(gca,'XLim',[20 140],'XDir','reverse','XTickLabel',[0,20,40,60,80,100,120],'YTick',[0.5,1],'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2 );
box off;
axis square;

plot(gca, [42.19+20, 42.19+20], [0, 1], 'k--', 'LineWidth',4);
plot(gca, [55.08+20, 55.08+20], [0, 1], 'k--','LineWidth',4);
plot(gca, [105.48+20, 105.48+20], [0, 1], 'k--','LineWidth',4);




