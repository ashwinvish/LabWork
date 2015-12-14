% TreeLength Plots

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
    %set(h, 'FaceColor', BarCMap);
    hold on;
    
end

set(gca,'XLim', [0 23] , 'XTick', 1:22, 'FontName', 'Arial', 'FontSize', 20);
set(gcf,'color','w');
xlabel('Neuron #', 'FontName', 'Arial', 'FontSize', 20);
ylabel('Number of branches', 'FontName', 'Arial', 'FontSize', 20);
box off;
axis square;
hold off;

figure();

MeanBranch = [ mean(BranchAlx), mean(BranchDbx), mean(BranchBarhl)];
SdBranch = [std(BranchAlx), std(BranchDbx), std(BranchBarhl)];

plot([ones(length(BranchAlx),1)' , 2*ones(length(BranchDbx),1)',3*ones(length(BranchBarhl),1)'] , [BranchAlx, BranchDbx,  BranchBarhl] ,'o', 'MarkerFaceColor', [0.7,0.7,0.7], 'MarkerEdgeColor','none', 'MarkerSize', 25 );
hold on;
plot(1,MeanBranch(1), 'o', 'MarkerFaceColor', calx, 'MarkerEdgeColor','none', 'MarkerSize', 35 );
plot(2,MeanBranch(2), 'o', 'MarkerFaceColor', cdbx, 'MarkerEdgeColor','none', 'MarkerSize', 35 );
plot(3,MeanBranch(3), 'o', 'MarkerFaceColor', cbarhl, 'MarkerEdgeColor','none', 'MarkerSize', 35 );
plot([1:3;1:3], [MeanBranch-SdBranch;MeanBranch+SdBranch], 'Color','k','LineWidth',2);

set(gca,'XTick', [1:3],'XTickLabel', {'group1'; 'group2'; 'group3'}, 'FontName', 'Arial', 'FontSize', 40);
set(gca, 'XLim', [0.5 3.5], 'FontName', 'Arial', 'FontSize', 40);
set(gcf,'color','w');
ylabel('Number of branches', 'FontName', 'Arial', 'FontSize', 40);
axis square;
box off;
hold off;

%% Plot length of trees

axLength = [];
clear temp;

for i = 1:size(cellIDs,2)
    if eval([cellIDs{i},'_axon'])>0
        for jj = 1:numel(eval([cellIDs{i},'_axon']))
            AxnNodes = eval([cellIDs{i},'_axon']);
            treeNodes(jj,1:3) = allTrees{i}{AxnNodes(jj)}{3};
        end
        lengthToPreNodeTest = findPathLength([cellIDs{i} , '_WithTags.swc'],[5,5,45],treeNodes);
        allLengthToPreNodeTest{i} = lengthToPreNodeTest;
        treeNodes = [];
        temp(1:size(allLengthToPreNodeTest{i},1)-1,i) = diff(sort(allLengthToPreNodeTest{i}));
        axLength = [axLength,sum(temp(:,i))];
        
    else
        axLength = [axLength,0];
        continue;
    end
end
denLength = cell2mat(allRawLength)-axLength;
sprintf('dendrite length / axon length = %d',sum(denLength)/sum(axLength));

%%

%figure();

[SortLength, I] = sort(cell2mat(allRawLength(:))/1000);
LengthRatio = denLength./axLength;
LengthRatio = LengthRatio(I);

calx = [1 0.5 0.3];
cdbx = [1 0.3 1];
cbarhl = [0.3 0.5 1];
LengthAlx = [];
LengthDbx = [];
LengthBarhl = [];
LengthRatioAlx = [];
LengthRatioDbx = [];
LengthRatioBarhl = [];

for i = 1:length(cellIDs)
    if ismember(cellIDs(I(i)),cellIDsAlx) == 1
        BarCMap = calx;
        LengthAlx = [LengthAlx,SortLength(i)];
        LengthRatioAlx = [LengthRatioAlx, LengthRatio(i)];
    elseif ismember(cellIDs(I(i)),cellIDsDbx) == 1
        BarCMap = cdbx;
        LengthDbx = [LengthDbx,SortLength(i)];
        LengthRatioDbx = [LengthRatioDbx, LengthRatio(i)];
    else
        BarCMap= cbarhl;
        LengthBarhl = [LengthBarhl,SortLength(i)];
        LengthRatioBarhl = [LengthRatioBarhl, LengthRatio(i)];
    end
    
    h1 = figure(1);
    ax1 = h1.CurrentAxes;
    bar(i,SortLength(i),'FaceColor', BarCMap);
    %set(h, 'FaceColor', BarCMap);
    hold on;
    
    h2 = figure(2);
    ax2 = h2.CurrentAxes
    plot(i, LengthRatio(i),'o','MarkerSize', 25,'MarkerFaceColor', BarCMap, 'MarkerEdgeColor','none');
    hold on;
    
end

set(ax1,'XLim', [0 23] , 'XTick', 1:5:22, 'FontName', 'Arial', 'FontSize', 20);
set(h1,'color','w');
xlabel(ax1,'Neuron #', 'FontName', 'Arial', 'FontSize', 20);
ylabel(ax1,'Total length in \mum', 'FontName', 'Arial', 'FontSize', 20);
box off;
axis square;
hold off;

set(ax2,'XLim', [0 23] , 'XTick', 1:5:22, 'FontName', 'Arial', 'FontSize', 20);
set(h2,'color','w');
xlabel(ax2,'Neuron #', 'FontName', 'Arial', 'FontSize', 20);
ylabel(ax2,'Ratio of dendrite to axon', 'FontName', 'Arial', 'FontSize', 20);
box off;
axis square;
hold off;




figure();
[ax,h1,h2] = plotyy(1:22,denLength(I)/1000,1:22,axLength(I)/1000, 'plot', 'plot');
h1.Marker =  'o';
h1.MarkerFaceColor = 'b';
h1.MarkerSize = 25;
h1.LineStyle = 'none';
%h1.MarkerEdgeColor = 'none';
h2.Marker = 'o';
h2.MarkerFaceColor =  'r';
h2.MarkerEdgeColor = 'none';
h2.MarkerSize = 25;
h2.LineStyle = 'none';
set(ax(1),'xcolor','k',  'FontName', 'Arial', 'FontSize', 20);
set(ax(2),'xcolor','k',  'FontName', 'Arial', 'FontSize', 20);
set(ax(1),'YLim', [0 500]);
set(ax(2),'YLim', [0 1100]);
set(ax(1),'ycolor','b');
set(ax(2),'ycolor','r');
xlabel('Neuron #',  'FontName', 'Arial', 'FontSize', 20);
ylabel(ax(1),'Dendritic length in \mum',  'FontName', 'Arial', 'FontSize', 20);
ylabel(ax(2),'Axonal length in \mum', 'FontName', 'Arial', 'FontSize', 20);
box off;
axis (ax(1), 'square');
axis (ax(2), 'square');
set(ax(1),'XLim',[1 23],'XTick', 1:5:22 );
set(ax(2),'XLim',[1 23],'XTick', 1:5:22);
set(gcf,'color','w');

figure();


