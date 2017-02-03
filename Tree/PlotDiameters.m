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
            AxnPlength = B;
            DenPlength = C;
        end
        AlxAxnDia(index,1:size(AxnDia{i},1)) = AxnDia{i}';
        AlxDenDia(index,1:size(DenDia{i},1)) = DenDia{i}';
        AlxAxonPlength(index,1:size(AxnPlength,1)) = AxnPlength;
        %AlxDenPlength(index,1:size(DenPlength,1)) = DenPlength;
        plot([1,3],[mean(AxnDia{i})/1000,mean(DenDia{i})/1000], '-o','MarkerFaceColor',[0.8,0.8,0.8],'MarkerSize', 35,'MarkerEdgeColor','k','Color','k', 'LineWidth',2);
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
            AxnPlength = B;
            DenPlength = C;
        end
        TransAxnDia(index,1:size(AxnDia{i},1)) = AxnDia{i}';
        TransDenDia(index,1:size(DenDia{i},1)) = DenDia{i}';
        TransAxonPlength(index,1:size(AxnPlength,1)) = AxnPlength;
       % TransDenPlength(index,1:size(DenPlength,1)) = DenPlength;
        plot([1,3],[mean(AxnDia{i})/1000,mean(DenDia{i})/1000], '-o','MarkerFaceColor',[0.8,0.8,0.8],'MarkerSize', 35,'MarkerEdgeColor','k','Color','k', 'LineWidth',2);
        hold all;
        index =index+1;
    end
end
pbaspect([1 2 1]);
plot([1,3],[mean(nonzeros(TransAxnDia))/1000,mean(nonzeros(TransDenDia))/1000],'-o','MarkerFaceColor','k','MarkerSize', 35,'MarkerEdgeColor','k','Color','k', 'LineWidth',2);

h = set(gca,'YLim',[0 0.35],'XLim',[0 4],'XTick',[1 3],'XTickLabel',{'axon', 'dendrite'},'Color',[ctrans 0.2], 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
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
            AxnPlength = B;
            DenPlength = C;
        end
        DbxAxnDia(index,1:size(AxnDia{i},1)) = AxnDia{i}';
        DbxDenDia(index,1:size(DenDia{i},1)) = DenDia{i}';
        DbxAxonPlength(index,1:size(AxnPlength,1)) = AxnPlength;
        %DbxDenPlength(index,1:size(DenPlength,1)) = DenPlength;
        plot([1,3],[mean(AxnDia{i})/1000,mean(DenDia{i})/1000], '-o','MarkerFaceColor',[0.8,0.8,0.8],'MarkerSize', 35,'MarkerEdgeColor','k','Color','k', 'LineWidth',2);
        hold all;
        index = index+1;
    end
    
end

pbaspect([1 2 1]);
plot([1,3],[mean(nonzeros(DbxAxnDia))/1000,mean(nonzeros(DbxDenDia))/1000],'-o','MarkerFaceColor','k','MarkerSize', 35,'MarkerEdgeColor','k','Color','k', 'LineWidth',2);

h = set(gca,'YLim',[0 0.35],'XLim',[0 4],'XTick',[1 3],'XTickLabel',{'axon', 'dendrite'},'Color',[cdbx 0.2], 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
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
            AxnPlength = B;
            DenPlength = C;
        end
        BarhlAxnDia(index,1:size(AxnDia{i},1)) = AxnDia{i}';
        BarhlDenDia(index,1:size(DenDia{i},1)) = DenDia{i}';
        BarhlAxonPlength(index,1:size(AxnPlength,1)) = AxnPlength;
        %BarhlDenPlength(index,1:size(DenPlength,1)) = DenPlength;
        plot([1,3],[nan,mean(DenDia{i})/1000], '-o','MarkerFaceColor',[0.8,0.8,0.8],'MarkerSize', 35,'MarkerEdgeColor','k','Color','k', 'LineWidth',2);
        hold all;
        index = index+1;
    end
    
end
pbaspect([1 2 1]);
plot([1,3],[nan,mean(nonzeros(BarhlDenDia))/1000],'-o','MarkerFaceColor','k','MarkerSize', 35,'MarkerEdgeColor','k','Color','k', 'LineWidth',2);

h = set(gca,'YLim',[0 0.35],'XLim',[0 4],'XTick',[1 3],'XTickLabel',{'axon', 'dendrite'},'Color',[cbarhl 0.2], 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
set(gcf,'color','w');
box off;
hold off;

%%

%alx taper
subplot(2,2,1);
GroupTaper(AlxAxonPlength,AlxDenPlength,AlxAxnDia,AlxDenDia);
subplot(2,2,2)
GroupTaper(TransAxonPlength,TransDenPlength,TransAxnDia,TransDenDia);
subplot(2,2,3)
GroupTaper(DbxAxonPlength,DbxDenPlength,DbxAxnDia,DbxDenDia);
subplot(2,2,4)
GroupTaper(BarhlAxonPlength, BarhlDenPlength, BarhlAxnDia, BarhlDenDia);
    
%% all diameters

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

% [h,p] = kstest2(nonzeros(AlxAxnDia),nonzeros(AlxDenDia));% if 1 then reject the null hypothesis that they are from the same distrubution with p
% [h,p] = kstest2(nonzeros(TransAxnDia),nonzeros(TransDenDia));
% [h,p] = kstest2(nonzeros(DbxAxnDia),nonzeros(DbxDenDia));
% [h,p] = kstest2(nonzeros(BarhlAxnDia),nonzeros(BarhlDenDia));

% Alx Axon and Den Dia stats
AllAlxAxnDia = [];
AllAlxDenDia = [];
for i = 1:size(AlxAxnDia,1)
    AlxAxnDiaMean(i) = mean(nonzeros(AlxAxnDia(i,:)));
    AlxAxnDiaSD(i) = std(nonzeros(AlxAxnDia(i,:)));
    AlxDenDiaMean(i) = mean(nonzeros(AlxDenDia(i,:)));
    AlxDenDiaSD(i) = std(nonzeros(AlxDenDia(i,:)))
    AllAlxAxnDia = [AllAlxAxnDia; nonzeros(AlxAxnDia(i,:))];
    AllAlxDenDia = [AllAlxDenDia; nonzeros(AlxDenDia(i,:))];
end

sprintf(' Alx axon diameter is %0.2f pm %0.2f', mean(AllAlxAxnDia)/1000, std(AllAlxAxnDia)/1000)
sprintf(' Alx den diameter is %0.2f pm %0.2f', mean(AllAlxDenDia)/1000, std(AllAlxDenDia)/1000)
[h,p] = ttest2( AllAlxAxnDia, AllAlxDenDia);
sprintf('Alx ttest %d, %f',h,p)

% Trans Axon and Den dia stats

AllTransAxnDia =[];
AllTransDenDia = [];

for i = 1:size(TransAxnDia,1)
    TransAxnDiaMean(i) = mean(nonzeros(TransAxnDia(i,:)));
    TransAxnDiaSD(i) = std(nonzeros(TransAxnDia(i,:)));
    TransDenDiaMean(i) = mean(nonzeros(TransDenDia(i,:)));
    TransDenDiaSD(i) = std(nonzeros(TransDenDia(i,:)));
    AllTransAxnDia = [AllTransAxnDia; nonzeros(TransAxnDia(i,:))];
    AllTransDenDia = [AllTransDenDia; nonzeros(TransDenDia(i,:))];
end

sprintf(' Trans axon diameter is %0.2f pm %0.2f', mean(AllTransAxnDia)/1000, std(AllTransAxnDia)/1000)
sprintf(' Trans den diameter is %0.2f pm %0.2f', mean(AllTransDenDia)/1000, std(AllTransDenDia)/1000)
[h,p] = ttest2(AllTransAxnDia, AllTransDenDia);
sprintf('Trans ttest %d, %f',h,p)

% Dbx axon and dendrite dia stats

AllDbxAxnDia = [];
AllDbxDenDia = [];

for i = 1:size(DbxAxnDia,1)
    DbxAxnDiaMean(i) = mean(nonzeros(DbxAxnDia(i,:)));
    DbxAxnDiaSD(i) = std(nonzeros(DbxAxnDia(i,:)));
    DbxDenDiaMean(i) = mean(nonzeros(DbxDenDia(i,:)));
    DbxDenDiaSD(i) = std(nonzeros(DbxDenDia(i,:)));
    AllDbxAxnDia =[AllDbxAxnDia; nonzeros(DbxAxnDia(i,:))];
    AllDbxDenDia = [AllDbxDenDia; nonzeros(DbxDenDia(i,:))];
end

sprintf(' Dbx axon diameter is %0.2f pm %0.2f', mean(AllDbxAxnDia)/1000, std(AllDbxAxnDia)/1000)
sprintf(' Dbx den diameter is %0.2f pm %0.2f', mean(AllDbxDenDia)/1000, std(AllDbxDenDia)/1000)
[h,p] = ttest2(AllDbxAxnDia, AllDbxDenDia);
sprintf('Dbx ttest %d, %f',h,p)

% Barhl axon and den stats

AllBarhlAxnDia = [];
AllBarhlDenDia = [];


for i = 1:size(BarhlAxnDia,1)
    BarhlAxnDiaMean(i) = mean(nonzeros(BarhlAxnDia(i,:)));
    BarhlAxnDiaSD(i) = std(nonzeros(BarhlAxnDia(i,:)))
    BarhlDenDiaMean(i) = mean(nonzeros(BarhlDenDia(i,:)));
    BarhlDenDiaSD(i) = std(nonzeros(BarhlDenDia(i,:)))
    AllBarhlAxnDia = [AllBarhlAxnDia ; nonzeros(BarhlAxnDia(i,:))];
    AllBarhlDenDia = [ AllBarhlDenDia; nonzeros(BarhlDenDia(i,:))];
end

sprintf(' Barhl axon diameter is %0.2f pm %0.2f', mean(AllBarhlAxnDia)/1000, std(AllBarhlAxnDia)/1000)
sprintf(' Barhl den diameter is %0.2f pm %0.2f', mean(AllBarhlDenDia)/1000, std(AllBarhlDenDia)/1000)
[h,p] = ttest2(AllBarhlAxnDia, AllBarhlDenDia);
sprintf('Barhl ttest %d, %f',h,p)

AllAxonsDia = [ AllAlxAxnDia; AllTransAxnDia; AllDbxAxnDia; AllBarhlAxnDia];
AllDendritesDia = [AllAlxDenDia; AllTransDenDia; AllDbxDenDia; AllBarhlDenDia];

AllAxnDiaMean = [AlxAxnDiaMean, TransAxnDiaMean, DbxAxnDiaMean, BarhlAxnDiaMean];
AllDenDiaMean = [AlxDenDiaMean, TransDenDiaMean, DbxDenDiaMean, BarhlDenDiaMean];
AllAxnSD = [AlxAxnDiaSD,TransAxnDiaSD, DbxAxnDiaSD, BarhlAxnDiaSD];
AllDenSd = [AlxDenDiaSD, TransDenDiaSD, DbxDenDiaSD, BarhlDenDiaSD];

figure();
subplot(2,1,1)
histogram(AllAxonsDia,'BinWidth',50,'FaceColor',[0,0.8,0]);
hold on;
plot([mean(AllAxonsDia) mean(AllAxonsDia)], ylim, 'k-','LineWidth',2);
set(gca,'XLim',[0 1400],'YLim', ylim, 'FontName', 'Arial', 'FontSize', 40);
subplot(2,1,2);
histogram(AllDendritesDia,'BinWidth',50,'FaceColor',[0.9,0,0]);
hold on;
plot([mean(AllDendritesDia) mean(AllDendritesDia)], ylim, 'k-','LineWidth',2);
xlabel('Diameter in nm', 'FontName', 'Arial', 'FontSize', 40);
set(gca,'XLim',[0 1400],'YLim', ylim,'FontName', 'Arial', 'FontSize', 40);



figure();
subplot(1,2,1);
plot(repmat([1;2],1,22), [AllAxnDiaMean; AllDenDiaMean], '-o', 'MarkerSize' , 25, 'MarkerFaceColor','auto', 'LineWidth', 2);
set(gca, 'XLim',[0.5 2.5], 'XTick',[1, 2], 'XTickLabel', {'Axon','Dendrite'}, 'FontName', 'Arial', 'FontSize', 40, 'LineWidth', 2);
ylabel('Diameter in nm', 'FontName', 'Arial', 'FontSize', 40);

subplot(1,2,2);
plot(repmat([1;2],1,2), [AllAxnDiaMean(find(AllAxnDiaMean>AllDenDiaMean)); AllDenDiaMean(find(AllDenDiaMean<AllAxnDiaMean))], '-o', 'MarkerSize' , 25, 'MarkerFaceColor','auto', 'LineWidth', 2);
set(gca, 'XLim',[0.5 2.5], 'YLim',[100 300],'XTick',[1, 2], 'XTickLabel', {'Axon','Dendrite'}, 'FontName', 'Arial', 'FontSize', 40, 'LineWidth', 2);
ylabel('Diameter in nm', 'FontName', 'Arial', 'FontSize', 40);



