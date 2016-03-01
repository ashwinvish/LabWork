% plot dendritic measures
figure();
index =1;
for i = 1:length(cellIDs)
    if ismember(cellIDs(i),cellIDsAlx) == 1
        [A,B,C,D,E] = CellDiameter(i,allTrees, cellIDs, false);
        if isempty(D)
            AxnDia{i} = [];
            DenDia{i} = E(:,4);
        else
            AxnDia{i} = D(:,4);
            DenDia{i} = E(:,4);
            
        end
        AlxAxnDia(index,1:size(AxnDia{i},1)) = AxnDia{i}';
        AlxDenDia(index,1:size(DenDia{i},1)) = DenDia{i}';
        plot([1,3],[mean(AxnDia{i})/1000,mean(DenDia{i})/1000], '-o','MarkerFaceColor',[0.8,0.8,0.8],'MarkerSize', 25,'MarkerEdgeColor','k','Color','k', 'LineWidth',2);
        hold all;
        index = index+1;
    end
end
plot([1,3],[mean(nonzeros(AlxAxnDia))/1000,mean(nonzeros(AlxDenDia))/1000],'-o','MarkerFaceColor','k','MarkerSize', 35,'MarkerEdgeColor','k','Color','k', 'LineWidth',2);

%plot([1;1],[mean(nonzeros(AlxAxnDia))/1000 + std(nonzeros(AlxAxnDia))/(sqrt(size(nonzeros(AlxAxnDia),1))*1000) ; mean(nonzeros(AlxAxnDia))/1000 - std(nonzeros(AlxAxnDia))/(sqrt(size(nonzeros(AlxAxnDia),1))*1000)], 'Color','k','LineWidth',2);
%plot([3;3],[mean(nonzeros(AlxDenDia))/1000 + std(nonzeros(AlxDenDia))/(sqrt(size(nonzeros(AlxDenDia),1))*1000); mean(nonzeros(AlxDenDia))/1000 - std(nonzeros(AlxDenDia))/(sqrt(size(nonzeros(AlxDenDia),1))*1000)], 'Color','k','LineWidth',2);
pbaspect([1 2 1]);
h = set(gca,'YLim',[0 0.35],'XLim',[0 4],'XTick',[1 3],'XTickLabel',{'axon', 'dendrite'},'Color',[calx 0.2], 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
set(gcf,'color','w');
ylabel('Diameter (\mum)', 'FontName', 'Arial', 'FontSize', 40);
box off;
hold off;


figure();
index = 1;
for i = 1:length(cellIDs)
    if ismember(cellIDs(i), cellIDsTrans) ==1
        [A,B,C,D,E] = CellDiameter(i,allTrees, cellIDs, false);
        if isempty(D)
            AxnDia{i} = [];
            DenDia{i} = E(:,4);
        else
            AxnDia{i} = D(:,4);
            DenDia{i} = E(:,4);
            
        end
        TransAxnDia(index,1:size(AxnDia{i},1)) = AxnDia{i}';
        TransDenDia(index,1:size(DenDia{i},1)) = DenDia{i}';
        plot([1,3],[mean(AxnDia{i})/1000,mean(DenDia{i})/1000], '-o','MarkerFaceColor',[0.8,0.8,0.8],'MarkerSize', 35,'MarkerEdgeColor','k','Color','k', 'LineWidth',2);
        hold all;
        index =index+1;
    end
end
pbaspect([1 2 1]);
plot([1,3],[mean(nonzeros(TransAxnDia))/1000,mean(nonzeros(TransDenDia))/1000],'-o','MarkerFaceColor','k','MarkerSize', 35,'MarkerEdgeColor','k','Color','k', 'LineWidth',2);

h = set(gca,'YLim',[0.1 0.35],'YTick',[],'YColor','none','XLim',[0 4],'XTick',[1 3],'XTickLabel',{'axon', 'dendrite'},'Color',[ctrans 0.2], 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
set(gcf,'color','w');
box off;
hold off;

figure();
index = 1;
for  i = 1:length(cellIDs)
    if ismember(cellIDs(i),cellIDsDbx) == 1
        [A,B,C,D,E] = CellDiameter(i,allTrees, cellIDs, false);
        if isempty(D)
            AxnDia{i} = [];
            DenDia{i} = E(:,4);
        else
            AxnDia{i} = D(:,4);
            DenDia{i} = E(:,4);
            
        end
        DbxAxnDia(index,1:size(AxnDia{i},1)) = AxnDia{i}';
        DbxDenDia(index,1:size(DenDia{i},1)) = DenDia{i}';
        plot([1,3],[mean(AxnDia{i})/1000,mean(DenDia{i})/1000], '-o','MarkerFaceColor',[0.8,0.8,0.8],'MarkerSize', 35,'MarkerEdgeColor','k','Color','k', 'LineWidth',2);
        hold all;
        index = index+1;
    end
    
end

pbaspect([1 2 1]);
plot([1,3],[mean(nonzeros(DbxAxnDia))/1000,mean(nonzeros(DbxDenDia))/1000],'-o','MarkerFaceColor','k','MarkerSize', 35,'MarkerEdgeColor','k','Color','k', 'LineWidth',2);

h = set(gca,'YLim',[0 0.35],'YTick',[],'YColor','none','XLim',[0 4],'XTick',[1 3],'XTickLabel',{'axon', 'dendrite'},'Color',[cdbx 0.2], 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
set(gcf,'color','w');
box off;
hold off;

figure();
index = 1;
for  i = 1:length(cellIDs)
    if ismember(cellIDs(i),cellIDsL) == 1
        [A,B,C,D,E] = CellDiameter(i,allTrees, cellIDs, false);
        
        if isempty(D)
            AxnDia{i} = [];
            DenDia{i} = E(:,4);
        else
            AxnDia{i} = D(:,4);
            DenDia{i} = E(:,4);
            
        end
        BarhlAxnDia(index,1:size(AxnDia{i},1)) = AxnDia{i}';
        BarhlDenDia(index,1:size(DenDia{i},1)) = DenDia{i}';
        plot([1,3],[nan,mean(DenDia{i})/1000], '-o','MarkerFaceColor',[0.8,0.8,0.8],'MarkerSize', 35,'MarkerEdgeColor','k','Color','k', 'LineWidth',2);
        hold all;
        index = index+1;
    end
    
end
pbaspect([1 2 1]);
plot([1,3],[nan,mean(nonzeros(BarhlDenDia))/1000],'-o','MarkerFaceColor','k','MarkerSize', 35,'MarkerEdgeColor','k','Color','k', 'LineWidth',2);

h = set(gca,'YLim',[0 0.35],'YTick',[],'YColor','none','XLim',[0 4],'XTick',[1 3],'XTickLabel',{'axon', 'dendrite'},'Color',[cbarhl 0.2], 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
set(gcf,'color','w');
box off;
hold off;


%% total number of 

edges = repmat({logspace(0,1,100)},1,22);

AxnHandles = cellfun(@histcounts,cellfun(@log,AxnDia, 'UniformOutput', false),edges, 'UniformOutput',false);
DenHandles = cellfun(@histcounts,cellfun(@log,DenDia, 'UniformOutput', false),edges, 'UniformOutput',false);

for i = 1:numel(cellIDs)
    AxnTotal(i,:) = AxnHandles{i};
    DenTotal(i,:) = DenHandles{i};
end

plot(logspace(0,1,100),[0,sum(AxnTotal)]./sum(sum(AxnTotal)),'-og', 'LineWidth',2);
hold on;
plot(logspace(0,1,100), [0,sum(DenTotal)]./sum(sum(DenTotal)),'-or', 'LineWidth',2);
ylabel('Proportion', 'FontName', 'Arial', 'FontSize', 40);
xlabel('Diameter in \mum', 'FontName', 'Arial', 'FontSize', 40);
axis square;
box off;
set(gca,'XScale','log','LineWidth',2,  'FontName', 'Arial', 'FontSize', 40);

%% Stats

% unique distributions; Kolmogorov-Simirinov test

[h,p] = kstest2(nonzeros(AlxAxnDia),nonzeros(AlxDenDia));% if 1 then reject the null hypothesis that they are from the same distrubution with p
[h,p] = kstest2(nonzeros(TransAxnDia),nonzeros(TransDenDia));
[h,p] = kstest2(nonzeros(DbxAxnDia),nonzeros(DbxDenDia));
[h,p] = kstest2(nonzeros(BarhlAxnDia),nonzeros(BarhlDenDia));

% Axon diameter is significantly different to Den diameter, ttest2

for i = 1:size(AlxAxnDia,1)
    AlxAxnDiaMean(i) = mean(nonzeros(AlxAxnDia(i,:)));
    AlxDenDiaMean(i) = mean(nonzeros(AlxDenDia(i,:)));
end

[h,p] = ttest2( AlxAxnDiaMean,AlxDenDiaMean);


for i = 1:size(TransAxnDia,1)
    TransAxnDiaMean(i) = mean(nonzeros(TransAxnDia(i,:)));
    TransDenDiaMean(i) = mean(nonzeros(TransDenDia(i,:)));
end

[h,p] = ttest2(TransAxnDiaMean,TransDenDiaMean);

for i = 1:size(DbxAxnDia,1)
    DbxAxnDiaMean(i) = mean(nonzeros(DbxAxnDia(i,:)));
    DbxDenDiaMean(i) = mean(nonzeros(DbxDenDia(i,:)));
end

[h,p] = ttest2(DbxAxnDiaMean,DbxDenDiaMean);

for i = 1:size(BarhlAxnDia,1)
    BarhlAxnDiaMean(i) = mean(nonzeros(BarhlAxnDia(i,:)));
    BarhlDenDiaMean(i) = mean(nonzeros(BarhlDenDia(i,:)));
end

[h,p] = ttest2(BarhlAxnDiaMean,BarhlDenDiaMean);
