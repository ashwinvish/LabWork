% TreeLength Plots Plot length of trees

%axLength = [];
%denLength = [];
clear temp;


figure();

index = 1;
        
for i = 1:numel(cellIDs)
    if ismember(cellIDs{i},cellIDsAlx) == 1
        a = axLength(i);
        d = denLength(i);
        AlxAxnLength(index) = a;
        AlxDenLenght(index) = d;
        plot([1,3],[a/1000,d/1000], '-o','MarkerFaceColor',[0.8,0.8,0.8],'MarkerSize', 35,'MarkerEdgeColor','k','Color','k', 'LineWidth',2);
        hold all;
        index = index+1;
    end
end
pbaspect([1 2 1]);
h = set(gca,'YLim',[-100 750],'XLim',[0 4],'XTick',[1 3],'XTickLabel',{'axon', 'dendrite'},'Color',[calx 0.2], 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
set(gcf,'color','w');
ylabel('Length in \mum', 'FontName', 'Arial', 'FontSize', 40);
box off;
hold off;
    
figure();
index = 1;
for i = 1:numel(cellIDs)
    if ismember(cellIDs{i},cellIDsTrans) == 1
        a = axLength(i);
        d = denLength(i);
        TransAxnLength(index) = a;
        TransDenLength(index) = d;
        index = index+1;
        plot([1,3],[a/1000,d/1000], '-o','MarkerFaceColor',[0.8,0.8,0.8],'MarkerSize', 35,'MarkerEdgeColor','k','Color','k', 'LineWidth',2);
        hold all;
    end
end
pbaspect([1 2 1]);
h = set(gca,'YLim',[-100 750],'YTick',[],'YColor','w','XLim',[0 4],'XTick',[1 3],'XTickLabel',{'axon', 'dendrite'},'Color',[ctrans 0.2], 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
set(gcf,'color','w');
ylabel('Length in \mum', 'FontName', 'Arial', 'FontSize', 40);
box off;
hold off;

figure();
index = 1;
for i = 1:numel(cellIDs)
    if ismember(cellIDs{i},cellIDsDbx) == 1
        a = axLength(i);
        d = denLength(i);
        DbxAxnLength(index) = a;
        DbxDenLength(index) = d;
        index = index+1;
        plot([1,3],[a/1000,d/1000], '-o','MarkerFaceColor',[0.8,0.8,0.8],'MarkerSize', 35,'MarkerEdgeColor','k','Color','k', 'LineWidth',2);
        hold all;
    end
end
pbaspect([1 2 1]);
h = set(gca,'YLim',[-100 750],'YTick',[],'YColor','w','XLim',[0 4],'XTick',[1 3],'XTickLabel',{'axon', 'dendrite'},'Color',[cdbx 0.2], 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
set(gcf,'color','w');
ylabel('Length in \mum', 'FontName', 'Arial', 'FontSize', 40);
box off;
hold off;

figure();
index = 1;
for i = 1:numel(cellIDs)
    if ismember(cellIDs{i},cellIDsL) == 1
        a = axLength(i);
        d = denLength(i);
        BarhlAxnLenght(index) = a;
        BarhlDenLength(index) = d;
        index = index+1;
        plot([1,3],[nan,d/1000], '-o','MarkerFaceColor',[0.8,0.8,0.8],'MarkerSize', 35,'MarkerEdgeColor','k','Color','k', 'LineWidth',2);
        hold all;
    end
end
pbaspect([1 2 1]);
h = set(gca,'YLim',[-100 750],'YTick',[],'YColor','w','XLim',[0 4],'XTick',[1 3],'XTickLabel',{'axon', 'dendrite'},'Color',[cbarhl 0.2], 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
set(gcf,'color','w');
ylabel('Length in \mum', 'FontName', 'Arial', 'FontSize', 40);
box off;
hold off;
    
%%

[SortLength, I] = sort(cell2mat(allRawLength(:))/1000);
LengthRatio = denLength./axLength;
LengthRatio = LengthRatio(I);

LengthAlx = [];
LengthTrans = [];
LengthDbx = [];
LengthBarhl = [];
LengthRatioAlx = [];
LengthRatioTrans = [];
LengthRatioDbx = [];
LengthRatioBarhl = [];

for i = 1:length(cellIDs)
    if ismember(cellIDs(I(i)),cellIDsAlx) == 1
        BarCMap = calx;
        LengthAlx = [LengthAlx,SortLength(i)];
        LengthRatioAlx = [LengthRatioAlx, LengthRatio(i)];
    elseif ismember(cellIDs(I(i)),cellIDsTrans) == 1
        BarCMap = ctrans
        LengthTrans = [LengthTrans,SortLength(i)];
        LengthRatioTrans = [LengthRatioTrans, LengthRatio(i)];
    elseif ismember(cellIDs(I(i)),cellIDsDbx) == 1
        BarCMap = cdbx;
        LengthDbx = [LengthDbx,SortLength(i)];
        LengthRatioDbx = [LengthRatioDbx, LengthRatio(i)];
    else
        BarCMap= cbarhl;
        LengthBarhl = [LengthBarhl,SortLength(i)];
        LengthRatioBarhl = [LengthRatioBarhl, LengthRatio(i)];
    end
      
end

figure();

MeanLength = [ mean(LengthAlx), mean(LengthTrans), mean(LengthDbx), mean(LengthBarhl)];
SdLength = [ std(LengthAlx), std(LengthTrans), std(LengthDbx), std(LengthBarhl)];

plot([ones(length(LengthAlx),1)' ,2*ones(length(LengthTrans),1)', 3*ones(length(LengthDbx),1)',4*ones(length(LengthBarhl),1)'] , [LengthAlx, LengthTrans ,LengthDbx ,LengthBarhl] ,'o', 'MarkerFaceColor', [0.7,0.7,0.7], 'MarkerEdgeColor','k', 'MarkerSize', 25 );
hold on;
plot(1,MeanLength(1), 'o', 'MarkerFaceColor', calx, 'MarkerEdgeColor','k', 'MarkerSize', 35 );
plot(2,MeanLength(2), 'o', 'MarkerFaceColor', ctrans, 'MarkerEdgeColor','k', 'MarkerSize', 35 );
plot(3,MeanLength(3), 'o', 'MarkerFaceColor', cdbx, 'MarkerEdgeColor','k', 'MarkerSize', 35 );
plot(4,MeanLength(4), 'o', 'MarkerFaceColor', cbarhl, 'MarkerEdgeColor','k', 'MarkerSize', 35 );
plot([1:4;1:4], [MeanLength-SdLength;MeanLength+SdLength], 'Color','k','LineWidth',2);

set(gca,'XTick', [1:4],'XTickLabel', {'group1'; 'group2'; 'group3'; 'group4'}, 'FontName', 'Arial', 'FontSize', 20);
set(gca, 'XLim', [0.5 4.5], 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
set(gcf,'color','w');
ylabel('Tree Length in \mum', 'FontName', 'Arial', 'FontSize', 40);
axis square;
box off;
hold off;




