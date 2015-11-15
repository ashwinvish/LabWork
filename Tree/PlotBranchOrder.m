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

calx = [1 0.5 0.3];
cdbx = [1 0.3 1];
cbarhl = [0.3 0.5 1];
BranchAlx = [];
BranchDbx = [];
BranchBarhl = [];

for i = 1:length(cellIDs)
    if ismember(cellIDs(I(i)),cellIDsAlx) == 1
        BarCMap = calx;
        BranchAlx = [BranchAlx,Branch(i)];
    elseif ismember(cellIDs(I(i)),cellIDsDbx) == 1
        BarCMap = cdbx;
        BranchDbx = [BranchDbx,Branch(i)];
    else
        BarCMap= cbarhl;
        BranchBarhl = [BranchBarhl,Branch(i)];
    end
    h = bar(i,Branch(i),'FaceColor', BarCMap);
    %[ax,h1,h2] = plotyy(i,Branch(i),i,allRawLength{I(i)}/1000,'bar', 'plot' );
    %set(h1, 'FaceColor', BarCMap);
    %set(h2,'Marker','*', 'color','k');
    hold on;
end

% set(gca,'XLim', [0 23] , 'XTick', 1:22, 'FontName', 'Arial', 'FontSize', 20);
% set(gcf,'color','w');
% xlabel('Neuron #', 'FontName', 'Arial', 'FontSize', 20);
% ylabel('Number of branches', 'FontName', 'Arial', 'FontSize', 20);
% box off;
% axis square;
% hold off;

figure();

MeanBranch = [ mean(BranchAlx), mean(BranchDbx), mean(BranchBarhl)];
SdBranch = [std(BranchAlx), std(BranchDbx), std(BranchBarhl)];

plot([ones(length(BranchAlx),1)' , 2*ones(length(BranchDbx),1)',3*ones(length(BranchBarhl),1)'] , [BranchAlx, BranchDbx,  BranchBarhl] ,'o', 'MarkerFaceColor', [0.7,0.7,0.7], 'MarkerEdgeColor','none', 'MarkerSize', 25 );
hold on;
plot(1,MeanBranch(1), 'o', 'MarkerFaceColor', calx, 'MarkerEdgeColor','none', 'MarkerSize', 25 );
plot(2,MeanBranch(2), 'o', 'MarkerFaceColor', cdbx, 'MarkerEdgeColor','none', 'MarkerSize', 25 );
plot(3,MeanBranch(3), 'o', 'MarkerFaceColor', cbarhl, 'MarkerEdgeColor','none', 'MarkerSize', 25 );
plot([1:3;1:3], [MeanBranch-SdBranch;MeanBranch+SdBranch], 'Color','k','LineWidth',2);

 set(gca,'XTick', [1:3],'XTickLabel', {'group1'; 'group2'; 'group3'}, 'FontName', 'Arial', 'FontSize', 20);
 set(gca, 'XLim', [0.5 3.5], 'FontName', 'Arial', 'FontSize', 20);
 set(gcf,'color','w');
 ylabel('Number of branches', 'FontName', 'Arial', 'FontSize', 20);
 axis square;
 box off;
hold off;



