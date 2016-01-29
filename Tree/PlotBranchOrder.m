% plot Branching order of all cells

figure();
h = tight_subplot(3,8,[.05 .05],[.05 .1],[.01 .01]);
for i =1:length(cellIDs)
    [Branches{i}, Terminals{i}, BranchOrder{i}] = TreeBranches(allTrees{i});
    axes(h(i));
    title(sprintf('Cell ID %s', cellIDs{i}));
    BranchOrderVisualizer(allTrees{i},[1],[BranchOrder{i}]);
end

figure();
h = tight_subplot(3,8,[.05 .05],[.05 .1],[.01 .01]);
for i =1:length(cellIDs)
    axes(h(i));
    histogram(BranchOrder{i},'BinLimits',[min(BranchOrder{i}), max(BranchOrder{i})]);
    title(sprintf('Cell ID %s', cellIDs{i}));
end

figure();

[rs,cs] = cellfun(@size,allSpine);
[Branch, I] = sort(cellfun(@sum,Branches)-rs);


BranchAlx = [];
BranchTrans = [];
BranchDbx = [];
BranchBarhl = [];

for i = 1:length(cellIDs)
    if ismember(cellIDs(I(i)),cellIDsAlx) == 1
        BarCMap = calx;
        BranchAlx = [BranchAlx,Branch(i)];
    elseif ismember(cellIDs(I(i)), cellIDsTrans) ==1
        BarCMap = ctrans;
        BranchTrans = [BranchTrans,Branch(i)];
    elseif ismember(cellIDs(I(i)),cellIDsDbx) == 1
        BarCMap = cdbx;
        BranchDbx = [BranchDbx,Branch(i)];
    else
        BarCMap= cbarhl;
        BranchBarhl = [BranchBarhl,Branch(i)];
    end
    plot(i,Branch(i),'o','MarkerFaceColor', BarCMap,'MarkerSize',35, 'MarkerEdgeColor','k');
    %[ax,h1,h2] = plotyy(i,Branch(i),i,allRawLength{I(i)}/1000,'bar', 'plot' );
    %set(h1, 'FaceColor', BarCMap);
    %set(h2,'Marker','*', 'color','k');
    hold on;
end

set(gca,'XLim', [0 23] , 'XTick', 1:22, 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
set(gcf,'color','w');
xlabel('Neuron #', 'FontName', 'Arial', 'FontSize', 40);
ylabel('Number of branches', 'FontName', 'Arial', 'FontSize', 40);
box off;
axis square;
hold off;

figure();

MeanBranch = [ mean(BranchAlx), mean(BranchTrans), mean(BranchDbx), mean(BranchBarhl)];
SdBranch = [std(BranchAlx),std(BranchTrans), std(BranchDbx), std(BranchBarhl)];

plot([ones(length(BranchAlx),1)' ,2*ones(length(BranchTrans),1)', 3*ones(length(BranchDbx),1)',4*ones(length(BranchBarhl),1)'] , [BranchAlx, BranchTrans, BranchDbx,  BranchBarhl] ,'o', 'MarkerFaceColor', [0.7,0.7,0.7], 'MarkerEdgeColor','k', 'MarkerSize', 25 );
hold on;
plot(1,MeanBranch(1), 'o', 'MarkerFaceColor', calx, 'MarkerEdgeColor','k', 'MarkerSize', 35 );
plot(2,MeanBranch(2), 'o', 'MarkerFaceColor', ctrans, 'MarkerEdgeColor','k', 'MarkerSize', 35 );
plot(3,MeanBranch(3), 'o', 'MarkerFaceColor', cdbx, 'MarkerEdgeColor','k', 'MarkerSize', 35 );
plot(4,MeanBranch(4), 'o', 'MarkerFaceColor', cbarhl, 'MarkerEdgeColor','k', 'MarkerSize', 35 );
plot([1:4;1:4], [MeanBranch-SdBranch;MeanBranch+SdBranch], 'Color','k','LineWidth',2);

set(gca,'XTick', [1:4],'XTickLabel', {'group1'; 'group2'; 'group3'; 'group4'}, 'FontName', 'Arial', 'FontSize', 20);
set(gca, 'XLim', [0.5 4.5], 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
set(gcf,'color','w');
ylabel('Number of branches', 'FontName', 'Arial', 'FontSize', 40);
axis square;
box off;
hold off;



