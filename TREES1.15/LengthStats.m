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

set(gca,'XLim', [0 23] , 'XTick', 1:22, 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
set(gcf,'color','w');
xlabel('Neuron #', 'FontName', 'Arial', 'FontSize', 40);
ylabel('Number of branches', 'FontName', 'Arial', 'FontSize', 40);
box off;
axis square;
hold off;

figure();

MeanBranch = [ mean(BranchAlx), mean(BranchDbx), mean(BranchBarhl)];
SdBranch = [std(BranchAlx), std(BranchDbx), std(BranchBarhl)];

plot([0.5*ones(length(BranchAlx),1)' , ones(length(BranchDbx),1)',1.5*ones(length(BranchBarhl),1)'] , [BranchAlx, BranchDbx,  BranchBarhl] ,'o', 'MarkerFaceColor', [0.7,0.7,0.7], 'MarkerEdgeColor','k', 'MarkerSize', 25 );
hold on;
plot(0.5,MeanBranch(1), 'o', 'MarkerFaceColor', calx, 'MarkerEdgeColor','k', 'MarkerSize', 35 );
plot(1,MeanBranch(2), 'o', 'MarkerFaceColor', cdbx, 'MarkerEdgeColor','k', 'MarkerSize', 35 );
plot(1.5,MeanBranch(3), 'o', 'MarkerFaceColor', cbarhl, 'MarkerEdgeColor','k', 'MarkerSize', 35 );
plot([0.5:0.5:1.5;0.5:0.5:1.5], [MeanBranch-SdBranch;MeanBranch+SdBranch], 'Color','k','LineWidth',2);

set(gca,'XTick', [0.5,1,1.5],'XTickLabel', {'group1'; 'group2'; 'group3'}, 'FontName', 'Arial', 'FontSize', 40,'LineWidth',2 );
set(gca, 'XLim', [0.25 1.75], 'FontName', 'Arial', 'FontSize', 40);
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
        AxnNodes = eval([cellIDs{i},'_axon']);
        tempLength = 0;
        for jj = 1:numel(eval([cellIDs{i},'_axon']))
            tempLength = tempLength + sum(allTrees{i}{AxnNodes(jj)}{1,4}{1,2});
        end
        axLength = [axLength,tempLength];
        
    else
        axLength = [axLength,0];
        continue;
    end
end
denLength = cell2mat(allRawLength)-axLength;
sprintf('dendrite length / axon length = %d',sum(denLength)/sum(axLength));


%%
% axLength = [];
% clear denLength;
% clear allLengthToPreNode;
% clear temp;
%
% for i = 1:size(cellIDs,2)
%     if eval([cellIDs{i},'_axon'])>0
%         for jj = 1:numel(eval([cellIDs{i},'_axon']))
%             AxnNodes = sort(eval([cellIDs{i},'_axon']));
%             treeNodes(jj,1:3) = allTrees{i}{AxnNodes(jj)}{1,3};
%         end
%         allLengthToPreNode{i} = findPathLength([cellIDs{i} , '_WithTags.swc'],allTrees{i},[5,5,45],treeNodes);
%         treeNodes = [];
%         %temp(1:size(allLengthToPreNode{i},1)-1,i) = diff(sort(allLengthToPreNode{i}));
%         %axLength = [axLength,sum(temp(:,i))];
%
%     else
%         %axLength = [axLength,0];
%         continue;
%     end
% end
% %denLength = cell2mat(allRawLength)-axLength;
% %sprintf('dendrite length / axon length = %d',sum(denLength)/sum(axLength));
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
hold on;

meanAlxLength = mean(LengthAlx);
meanDbxLength = mean(LengthDbx);
meanBarhlLength = mean(LengthBarhl);

stdAxlLength = std(LengthAlx);
stdDbxLength = std(LengthDbx);
stdBarhlLength = std(LengthBarhl);

meanLength = [meanAlxLength,meanDbxLength,meanBarhlLength ];
stdLength = [stdAxlLength,stdDbxLength,stdBarhlLength];




plot(0.5*ones(size(LengthAlx,2)), LengthAlx,'o', 'MarkerFaceColor', [0.7,0.7,0.7], 'MarkerEdgeColor','k', 'MarkerSize', 25);
plot(ones(size(LengthDbx,2)), LengthDbx,'o', 'MarkerFaceColor', [0.7,0.7,0.7], 'MarkerEdgeColor','k', 'MarkerSize', 25);
plot(1.5*ones(size(LengthBarhl,2)), LengthBarhl,'o', 'MarkerFaceColor', [0.7,0.7,0.7], 'MarkerEdgeColor','k', 'MarkerSize', 25);
plot(0.5,meanAlxLength, 'o', 'MarkerFaceColor', calx, 'MarkerEdgeColor','k', 'MarkerSize', 35 );
plot(1,meanDbxLength, 'o', 'MarkerFaceColor', cdbx, 'MarkerEdgeColor','k', 'MarkerSize', 35 );
plot(1.5,meanBarhlLength, 'o', 'MarkerFaceColor', cbarhl, 'MarkerEdgeColor','k', 'MarkerSize', 35 );
plot([0.5:0.5:1.5;0.5:0.5:1.5], [meanLength-stdLength;meanLength+stdLength], 'Color','k','LineWidth',2);

set(gca,'XTick', [0.5,1,1.5],'XTickLabel', {'group1'; 'group2'; 'group3'}, 'FontName', 'Arial', 'FontSize', 40,'LineWidth',2 );
set(gca, 'XLim', [0.25 1.75], 'FontName', 'Arial', 'FontSize', 40);
set(gcf,'color','w');
ylabel('Tree Length in \mum', 'FontName', 'Arial', 'FontSize', 40);
axis square;
box off;
hold off;



figure();
% [ax,h1,h2] = plotyy(1:22,denLength(I)/1000,1:22,axLength(I)/1000, 'plot', 'plot');
% h1.Marker =  'o';
% h1.MarkerFaceColor = 'b';
% h1.MarkerSize = 25;
% h1.LineStyle = 'none';
% %h1.MarkerEdgeColor = 'none';
% h2.Marker = 'o';
% h2.MarkerFaceColor =  'r';
% h2.MarkerEdgeColor = 'none';
% h2.MarkerSize = 25;
% h2.LineStyle = 'none';
% set(ax(1),'xcolor','k',  'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
% set(ax(2),'xcolor','k',  'FontName', 'Arial', 'FontSize', 40,'LineWidth',2);
% set(ax(1),'YLim', [0 500], 'YTick', [0,125,250,375,500]);
% set(ax(2),'YLim', [0 1100], 'YTick',[0,250,500, 750,1100]);
% set(ax(1),'ycolor','b');
% set(ax(2),'ycolor','r');
% xlabel('Neuron #',  'FontName', 'Arial', 'FontSize', 40);
% ylabel(ax(1),'Dendritic length in \mum',  'FontName', 'Arial', 'FontSize', 40);
% ylabel(ax(2),'Axonal length in \mum', 'FontName', 'Arial', 'FontSize', 40);
% box off;
% axis (ax(1), 'square');
% axis (ax(2), 'square');
% set(ax(1),'XLim',[1 23],'XTick', 1:5:22 );
% set(ax(2),'XLim',[1 23],'XTick', 1:5:22);
% set(gcf,'color','w');


h1 = plot(1:22, denLength(I)/1000 );
h1.Marker =  'o';
h1.MarkerFaceColor = 'b';
h1.MarkerSize = 25;
h1.LineStyle = 'none';
hold on;
h2 = plot(1:22 , axLength(I)/1000);
h2.Marker = 'o';
h2.MarkerFaceColor =  'r';
h2.MarkerEdgeColor = 'none';
h2.MarkerSize = 25;
h2.LineStyle = 'none';
xlabel('Neuron #',  'FontName', 'Arial', 'FontSize', 40);
ylabel('Length in \mum' , 'FontName', 'Arial', 'FontSize', 40);
set(gca,'XLim',[0,22],'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
set(gcf,'color','w');
axis square;
box off;

%%
figure();
for i = 1:1:6
    plot([1,2], [denLength(I(i))/1000, axLength(I(i))], '-ko','MarkerSize',25, 'MarkerFaceColor', [0.7,0.7,0.7], 'LineWidth',2);
    hold on;
end
plot([1,2], [mean(denLength(I(1:6)))/1000, mean(axLength(I(1:6)))/1000], '-ko','MarkerSize',35,'MarkerFaceColor', 'k');
set(gca,'YLim', [0, 800], 'XLim', [0.5,2.5], 'XTick', [1,2],'XTickLabel',{'Dendrite','Axon'},'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
box off;

figure();
for  j = 7:16
    plot([1,2], [denLength(I(j))/1000, axLength(I(j))/1000], '-ko','MarkerSize',25, 'MarkerFaceColor', [0.7,0.7,0.7], 'LineWidth',2);
    hold on;
end
plot([1,2], [mean(denLength(I(7:13)))/1000, mean(axLength(I(7:13)))/1000], '-ko','MarkerSize',35,'MarkerFaceColor', 'k');
set(gca,'YLim', [0, 800], 'XLim', [0.5,2.5], 'XTick', [1,2],'XTickLabel',{'Dendrite','Axon'},'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2 );
box off;  

figure();
for  k = 17:22
    plot([1,2], [denLength(I(k))/1000, axLength(I(k))/1000], '-ko','MarkerSize',25, 'MarkerFaceColor', [0.7,0.7,0.7], 'LineWidth',2);
    hold on;
end
plot([1,2], [mean(denLength(I(14:22)))/1000, mean(axLength(I(14:22)))/1000], '-ko','MarkerSize',35,'MarkerFaceColor', 'k');
set(gca,'YLim', [0, 800], 'XLim', [0.5,2.5], 'XTick', [1,2],'XTickLabel',{'Dendrite','Axon'},'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2 );
box off;



