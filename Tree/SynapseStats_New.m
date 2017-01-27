% Synapse Stats
AlxPre = [];
AlxPreLength = [];
AlxPost = [];
AlxPostLength = [];
TransPre = [];
TransPreLength = [];
TransPost = [];
TransPostLength = [];
DbxPost = [];
DbxPostLength = [];
BarhlPost = [];
BarhlPostLength = [];

for kk = 1:numel(cellIDs)
    % subplot(3,8,kk);
    if ismember(cellIDs{kk},cellIDsAlx)==1
        AlxPost  = [AlxPost; allPost{kk}];
        AlxPre = [AlxPre; allPreSynapse{kk}];
        AlxPostLength = [AlxPostLength, length(allPost{kk})];
        AlxPreLength = [AlxPreLength, length(allPreSynapse{kk})];
        
    elseif ismember(cellIDs{kk},cellIDsTrans)==1
        TransPost  = [TransPost; allPost{kk}];
        TransPre = [TransPre; allPreSynapse{kk}];
        TransPostLength = [TransPostLength, length(allPost{kk})];
        TransPreLength = [TransPreLength, length(allPreSynapse{kk})];
        
    elseif ismember(cellIDs{kk},cellIDsDbx)==1
        DbxPost = [DbxPost; allPost{kk}];
        DbxPostLength = [DbxPostLength, length(allPost{kk})];
        
    else
        BarhlPost = [BarhlPost; allPost{kk}];
        BarhlPostLength = [BarhlPostLength, length(allPost{kk})];
    end
end

% bar plot of number of Postsynapses
figure();
MeanPost = [mean(AlxPostLength), mean(TransPostLength), mean(DbxPostLength), mean(BarhlPostLength)];
StdPost = [std(AlxPostLength), std(TransPostLength), std(DbxPostLength), std(BarhlPostLength)];
MeanPre = [mean(AlxPreLength), mean(TransPreLength)];
StdPre = [std(AlxPreLength), std(TransPreLength)];


plot([ones(length(AlxPostLength),1)',2*ones(length(TransPostLength),1)', 3*ones(length(DbxPostLength),1)',4*ones(length(BarhlPostLength),1)'],[AlxPostLength,TransPostLength,DbxPostLength,BarhlPostLength], 'Marker', 'o', 'MarkerFaceColor', [0.7,0,0],'MarkerSize',35,'LineStyle','none', 'MarkerEdgeColor','w');
hold on;
plot(1.2*ones(length(AlxPreLength),1)',AlxPreLength, 'Marker', 'o', 'MarkerFaceColor', [0,0.7,0],'MarkerSize',35,'LineStyle','none',  'MarkerEdgeColor','w');
plot(2.2*ones(length(TransPreLength),1)',TransPreLength, 'Marker', 'o', 'MarkerFaceColor', [0,0.7,0],'MarkerSize',35,'LineStyle','none',  'MarkerEdgeColor','w' );

plot(1.2,MeanPre(1),'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor','w', 'MarkerSize', 35 );
plot(2.2,MeanPre(2),'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor','w', 'MarkerSize', 35 );

plot([1.2,2.2;1.2,2.2], [MeanPre-StdPre;MeanPre+StdPre], 'Color','k','LineWidth',2);
plot(1,MeanPost(1), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor','w', 'MarkerSize', 35 );
plot(2,MeanPost(2), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor','w', 'MarkerSize', 35 );
plot(3,MeanPost(3), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor','w', 'MarkerSize', 35 );
plot(4,MeanPost(4), 'o', 'MarkerFaceColor',' k', 'MarkerEdgeColor','w', 'MarkerSize', 35 );
plot([1:4;1:4], [MeanPost-StdPost;MeanPost+StdPost], 'Color','k','LineWidth',2);

set(gca, 'XTick', [1:4],'XTickLabel', {'group1'; 'group2'; 'group3'; 'group4'},'XLim', [0.5 4.5], 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
ylabel('Number of sites', 'FontName', 'Arial', 'FontSize', 40);
set(gcf,'color','w');
legend('Postsynapse','Presynapse');
legend('boxoff');
box off;
axis square;

cpost = [0.8,0,0];
cpre = [0,0.7,0];

figure(2);

h = boxplot([AlxPostLength'; AlxPreLength'; TransPostLength';TransPreLength'; DbxPostLength'; BarhlPostLength'],...
    [ones(size(AlxPostLength,2),1); 2*ones(size(AlxPreLength,2),1); 3*ones(size(TransPostLength,2),1); 4*ones(size(TransPreLength,2),1); ...
    5*ones(size(DbxPostLength,2),1); 6*ones(size(BarhlPostLength,2),1)],'Notch','off', 'Symbol', 'ko',...
    'Colors',[cpost;cpre;cpost;cpre;cpost;cpost],...
    'OutlierSize', 25);
set(h,{'linew'},{4});
h1 = findobj(gca,'tag','Median');
set(h1,'Color','k');
set(gca,'FontName', 'Arial', 'FontSize', 40, 'LineWidth',4);
box off;

hold on;

plot([ones(size(AlxPostLength,2),1)', 2*ones(size(AlxPreLength,2),1)', 3*ones(size(TransPostLength,2),1)', 4*ones(size(TransPreLength,2),1)', ...
    5*ones(size(DbxPostLength,2),1)', 6*ones(size(BarhlPostLength,2),1)'],[AlxPostLength, AlxPreLength, TransPostLength,TransPreLength, DbxPostLength, BarhlPostLength], ...
     'Marker', 'o', 'MarkerFaceColor','none' ,'MarkerSize',25,'LineStyle','none', 'MarkerEdgeColor','k', 'LineWidth', 4);
axis square;


[p,h] = ranksum(AlxPostLength,DbxPostLength);



%% Distribution of inter-synaptic distance

clear InterPostSynapticDistance;
InterPostSynapticDistance = [];
figure;

for ii = 1:numel(cellIDs)
    InterPostSynapticDistance = [InterPostSynapticDistance;nonzeros(allPostDiff{ii})];
end
subplot(2,1,1)
histogram(InterPostSynapticDistance/1000, 'FaceColor',[0.9,0,0],'BinWidth',500/1000, 'Normalization','probability');
hold on;
set(gca,'XTick',[0,5,10,15,20,25],'XLim',[-0.5,25], 'YLim', [-0.1 0.6],'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);


box off;
set(gcf,'color','w');

clear x_val;
clear y_val;

clear InterPreSynapticDistance;
InterPreSynapticDistance = [];
tempPre =[];

for ii = 1:size(cellIDs,2)
    if ~isempty(allPreDiff{ii})
        InterPreSynapticDistance = [InterPreSynapticDistance; nonzeros(allPreDiff{ii})];
    end
end
subplot(2,1,2)
histogram(InterPreSynapticDistance/1000, 'FaceColor',[0,0.8,0],'BinWidth',0.5, 'Normalization','Probability');
hold on;
xlabel('Inter synaptic distance (\mum)', 'FontName', 'Arial', 'FontSize', 40);
ylabel('Probability', 'FontName', 'Arial', 'FontSize', 40);
set(gca, 'XTick',[0,5,10,15,20,25],'XLim',[-0.5,25], 'YLim', [-0.1 0.6],'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
box off;
set(gcf,'color','w');

figure();

h1 = cdfplot(InterPostSynapticDistance/1000 );
h1.Color = [0.9,0,0];
h1.LineWidth = 2;
hold on;
h2 = cdfplot(InterPreSynapticDistance/1000);
h2.Color = [0,0.8,0];
h2.LineWidth = 2;
legend({'Den.', 'Axon'},'FontName', 'Arial', 'FontSize', 40, 'Location','southeast' )
xlabel('Inter synaptic distance (\mum)', 'FontName', 'Arial', 'FontSize', 40);
ylabel('F(x)', 'FontName', 'Arial', 'FontSize', 40);
set(gca,'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
set(gcf,'color','w');
box off;
axis square


figure();
histogram(InterPostSynapticDistance/1000, 'FaceColor',[0.9,0.0,0],'BinWidth',0.5, 'Normalization','probability');
hold on;
line([mean(InterPostSynapticDistance/1000),mean(InterPostSynapticDistance/1000)], [0, 0.6], 'color', [0.9,0,0], 'LineWidth', 2)
histogram(InterPreSynapticDistance/1000, 'FaceColor',[0,0.8,0],'BinWidth',0.5, 'Normalization','probability');
line([mean(InterPreSynapticDistance/1000),mean(InterPreSynapticDistance/1000)], [0, 0.6], 'color', [0,0.8,0], 'LineWidth', 2)

axis square;
xlabel('Inter synaptic distance (\mum)', 'FontName', 'Arial', 'FontSize', 40);
ylabel('Probability', 'FontName', 'Arial', 'FontSize', 40);
set(gca, 'XTick',[0,5,10,15,20,25],'XLim',[-0.5,25],'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
set(gcf,'color','w');
box off;
axis square
%%
AlxInterPost = [];
TransInterPost = [];
DbxInterPost = [];
BarhlInterPost = [];


for i = 1:numel(cellIDs)
    if ismember(cellIDs{i}, cellIDsAlx) ==1
        AlxInterPost = [AlxInterPost; nonzeros(allPostDiff{i})];
    elseif ismember(cellIDs{i}, cellIDsTrans) ==1
        TransInterPost = [TransInterPost; nonzeros(allPostDiff{i})];
    elseif ismember(cellIDs{i}, cellIDsDbx) ==1
        DbxInterPost = [DbxInterPost; nonzeros(allPostDiff{i})];
    else
        BarhlInterPost = [BarhlInterPost; nonzeros(allPostDiff{i})];
    end
end




AlxInterPre = [];
TransInterPre = [];
DbxInterPre = [];
BarhlInterPre = [];


for i = 1:numel(cellIDs)
    if ismember(cellIDs{i}, cellIDsAlx) ==1
        AlxInterPre = [AlxInterPre; nonzeros(allPreDiff{i})];
    elseif ismember(cellIDs{i}, cellIDsTrans) ==1
        TransInterPre = [TransInterPre; nonzeros(allPreDiff{i})];
    elseif ismember(cellIDs{i}, cellIDsDbx) ==1
        DbxInterPre = [DbxInterPre; nonzeros(allPreDiff{i})];
    else
        BarhlInterPre = [BarhlInterPre; nonzeros(allPreDiff{i})];
    end
end

figure(1);
subplot(2,4,1);
histogram(AlxInterPost/1000, 'FaceColor',[0.9,0,0], 'BinWidth',1, 'Normalization','probability', 'EdgeColor', 'none');
set(gca, 'XLim',[-2, 40], 'YLim',[-0.05,0.6],'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
box off;
set(gcf,'color','w');
set(gca,'color',[calx,0.2]);

subplot(2,4,5);
histogram(AlxInterPre/1000, 'FaceColor',[0,0.8,0],'BinWidth',1,'Normalization','probability', 'EdgeColor', 'none');
%ylabel('Probability', 'FontName', 'Arial', 'FontSize', 40);
set(gca, 'XLim',[-2, 40], 'YLim',[-0.05,0.6],'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
box off;
set(gcf,'color','w');
set(gca,'color',[calx,0.2]);
%axis square;


subplot(2,4,2);
histogram(TransInterPost/1000, 'FaceColor',[0.9,0,0], 'BinWidth',1, 'Normalization','probability', 'EdgeColor', 'none');
set(gca, 'XLim',[-2, 40], 'YLim',[-0.05,0.6],'YColor','k','FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
box off;
set(gcf,'color','w');
set(gca,'color',[ctrans,0.2]);

subplot(2,4,6);
histogram(TransInterPre/1000,'FaceColor',[0,0.8,0],'BinWidth',1,'Normalization','probability', 'EdgeColor', 'none');
ylabel('Probability', 'FontName', 'Arial', 'FontSize', 40);
set(gca, 'XLim',[-2, 40], 'YLim',[-0.05,0.6],'YColor','k','FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
box off;
set(gcf,'color','w');
set(gca,'color',[ctrans,0.2]);
%axis square;

subplot(2,4,3);
histogram(DbxInterPost/1000,  'FaceColor',[0.9,0,0], 'BinWidth',1, 'Normalization','probability', 'EdgeColor', 'none');
%ylabel('Probability', 'FontName', 'Arial', 'FontSize', 40);
set(gca, 'XLim',[-2, 40], 'YLim',[-0.05,0.6],'YColor','k','FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
box off;
set(gcf,'color','w');
set(gca,'color',[cdbx,0.2]);
%axis square;

subplot(2,4,4);
histogram(BarhlInterPost/1000, 'FaceColor',[0.9,0,0], 'BinWidth',1, 'Normalization','probability', 'EdgeColor', 'none');
xlabel('Synaptic pathlenght (\mum)', 'FontName', 'Arial', 'FontSize', 40);
%ylabel('Probability', 'FontName', 'Arial', 'FontSize', 40);
set(gca, 'XLim',[-2, 40], 'YLim',[-0.05,0.6],'YColor','k','FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
box off;
set(gcf,'color','w');
set(gca,'color',[cbarhl,0.2]);
%axis square;




pm = char(177);
figure(2);

subplot(2,2,1);
h1 = cdfplot(AlxInterPre/1000);
h1.Color = [0,0.8,0];
h1.LineWidth = 4;
hold on;

h1 = cdfplot(AlxInterPost/1000);
h1.Color = [0.9,0,0];
h1.LineWidth = 4;

text(20,0.8,sprintf('%1.2f%c%1.2f', mean(AlxInterPost/1000),pm, std(AlxInterPost/1000)), 'FontName', 'Arial', 'FontSize', 40, 'color', [0.9,0,0]);
text(20,0.7,sprintf('CV = %1.2f',std(AlxInterPost)/mean(AlxInterPost)), 'FontName', 'Arial', 'FontSize', 40, 'color', [0.9,0,0]);
text(20,0.6,sprintf('%1.2f%c%1.2f', mean(AlxInterPre/1000),pm, std(AlxInterPre/1000)), 'FontName', 'Arial', 'FontSize', 40, 'color', [0,0.8,0]);
text(20,0.5, sprintf('CV = %1.2f',std(AlxInterPre)/mean(AlxInterPre)), 'FontName', 'Arial', 'FontSize', 40, 'color', [0,0.8,0]);

xlabel('Pathlength (\mum)','FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
ylabel('Cumulative fraction', 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
set(gca,'Color', [calx,0.2], 'XLim',[0, 40],'YTick',[0,0.2,0.4,0.6,0.8,1],'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
axis square;
box off;



subplot(2,2,2);

h2 = cdfplot(TransInterPre/1000);
h2.Color = [0,0.8,0];
h2.LineWidth = 4;
hold on;

h2 = cdfplot(TransInterPost/1000);
h2.Color = [0.9,0,0];
h2.LineWidth = 4;

text(20,0.8,sprintf('%1.2f%c%1.2f', mean(TransInterPost/1000),pm, std(TransInterPost/1000)), 'FontName', 'Arial', 'FontSize', 40, 'color', [0.9,0,0]);
text(20,0.7,sprintf('CV = %1.2f',std(TransInterPost)/mean(TransInterPost)), 'FontName', 'Arial', 'FontSize', 40, 'color', [0.9,0,0]);
text(20,0.6,sprintf('%1.2f%c%1.2f', mean(TransInterPre/1000),pm, std(TransInterPre/1000)), 'FontName', 'Arial', 'FontSize', 40, 'color', [0,0.8,0]);
text(20,0.5, sprintf('CV = %1.2f',std(TransInterPre)/mean(TransInterPre)), 'FontName', 'Arial', 'FontSize', 40, 'color', [0,0.8,0]);


xlabel('Pathlength (\mum)','FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
ylabel('Cumulative fraction', 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
set(gca,'Color', [ctrans,0.2], 'XLim',[0, 40],'YTick',[0,0.2,0.4,0.6,0.8,1],'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
axis square;
box off;


subplot(2,2,3);

h3 = cdfplot(DbxInterPost/1000);
h3.Color = [0.9,0,0];
h3.LineWidth = 4;

text(20,0.8,sprintf('%1.2f%c%1.2f', mean(DbxInterPost/1000),pm, std(DbxInterPost/1000)), 'FontName', 'Arial', 'FontSize', 40, 'color', [0.9,0,0]);
text(20,0.7,sprintf('CV = %1.2f',std(DbxInterPost)/mean(DbxInterPost)), 'FontName', 'Arial', 'FontSize', 40, 'color', [0.9,0,0]);

xlabel('Pathlength (\mum)','FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
ylabel('Cumulative fraction', 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
set(gca,'Color', [cdbx, 0.2], 'XLim',[0, 40],'YTick',[0,0.2,0.4,0.6,0.8,1],'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
axis square;
box off;


subplot(2,2,4);
h4 = cdfplot(BarhlInterPost/1000);
h4.Color = [0.9,0,0];
h4.LineWidth = 4;

text(20,0.8,sprintf('%1.2f%c%1.2f', mean(BarhlInterPost/1000),pm, std(BarhlInterPost/1000)), 'FontName', 'Arial', 'FontSize', 40, 'color', [0.9,0,0]);
text(20,0.7,sprintf('CV = %1.2f',std(BarhlInterPost)/mean(BarhlInterPost)), 'FontName', 'Arial', 'FontSize', 40, 'color', [0.9,0,0]);

xlabel('Pathlength (\mum)','FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
ylabel('Cumulative fraction', 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
set(gca,'Color', [cbarhl, 0.2], 'XLim',[0, 40],'YTick',[0,0.2,0.4,0.6,0.8,1],'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
axis square;
box off;

figure(3);

h1 = cdfplot(AlxInterPost/1000);
h1.Color = calx;
h1.LineWidth = 4;

hold on;

h2 = cdfplot(TransInterPost/1000);
h2.Color = ctrans;
h2.LineWidth = 4;


h3 = cdfplot(DbxInterPost/1000);
h3.Color = cdbx;
h3.LineWidth = 4;

h4 = cdfplot(BarhlInterPost/1000);
h4.Color = cbarhl;
h4.LineWidth = 4;

axis square


figure(4)

h1 = cdfplot(AlxInterPre/1000);
h1.Color = calx;
h1.LineWidth = 4;

hold on;

h2 =  cdfplot(TransInterPre/1000);
h2.Color = ctrans;
h2.LineWidth = 4;


figure(5);

h = boxplot([AlxInterPost/1000; TransInterPost/1000; DbxInterPost/1000; BarhlInterPost/1000]...
, [ones(size(AlxInterPost,1),1); 2*ones(size(TransInterPost,1),1); 3*ones(size(DbxInterPost,1),1); 4*ones(size(BarhlInterPost,1),1)],...
'Notch','off', 'Symbol', 'ko', 'Labels',{'Ipsi', 'Ipsi-Contra','Contra', 'Unknown'},'Colors',[calx;ctrans;cdbx; cbarhl] ...
,'OutlierSize', 25);
set(h,{'linew'},{4});
h1 = findobj(gca,'tag','Median');
set(h1,'Color','k');
set(gca,'YLim', [0,10],'FontName', 'Arial', 'FontSize', 40, 'LineWidth',4);
box off;

% hold on;
% plot([ones(size(AlxInterPost,1),1)', 2*ones(size(TransInterPost,1),1)', 3*ones(size(DbxInterPost,1),1)', 4*ones(size(BarhlInterPost,1),1)'],...
%    [AlxInterPost/1000; TransInterPost/1000; DbxInterPost/1000; BarhlInterPost/1000],'o','MarkerSize', 25, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'k', 'LineWidth',2 );

hold off;
axis square;

InterPrePValues = []



figure(6);

h = boxplot([AlxInterPre/1000; TransInterPre/1000],[ones(size(AlxInterPre,1),1); 2*ones(size(TransInterPre,1),1)],...
'Notch','off', 'Symbol', 'ko', 'Labels',{'Ipsi', 'Ipsi-Contra'},'Colors',[calx;ctrans],'OutlierSize', 25);
set(h,{'linew'},{4});
h1 = findobj(gca,'tag','Median');
set(h1,'Color','k');
set(gca,'YLim', [0,10],'FontName', 'Arial', 'FontSize', 40, 'LineWidth',4);

box off;




%% Distribution of pre and post synaptic path lenghts

allPostSynapticLength = [];
allPreSynapticLength = [];


for i = 1:size(cellIDs,2)
    allPostSynapticLength = [allPostSynapticLength;allLengthToPostNode{i}];
end

for i = 1:size(cellIDs,2)
    allPreSynapticLength = [allPreSynapticLength;allLengthToPreNode{i}];
end

figure();
histogram(allPostSynapticLength/1000,'FaceColor',[0.9,0,0],'BinWidth',10,'Normalization','probability'); % dimensions in microns
hold on;
line([mean(allPostSynapticLength/1000),mean(allPostSynapticLength/1000)], [0,0.1], 'Color',[0.9,0,0], 'LineWidth', 4);
text(mean(allPostSynapticLength/1000)+10, 0.1, num2str(mean(allPostSynapticLength/1000)),'FontName', 'Arial', 'FontSize', 40 );
histogram(allPreSynapticLength/1000,'FaceColor',[0,0.8,0],'BinWidth',10,'Normalization','probability'); % dimensions in microns
line([mean(allPreSynapticLength/1000),mean(allPreSynapticLength/1000)], [0,0.1], 'Color',[0,0.8,0], 'LineWidth', 4);
text(mean(allPreSynapticLength/1000)+10, 0.1, num2str(mean(allPreSynapticLength/1000)),'FontName', 'Arial', 'FontSize', 40 );

%legend({'PostSynapse','PreSynapse'}, 'box','off');
xlabel('Synaptic pathlenght (\mum)', 'FontName', 'Arial', 'FontSize', 40);
ylabel('Count', 'FontName', 'Arial', 'FontSize', 40);
set(gca,'XLim',[0,350], 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
box off;
set(gcf,'color','w');
axis square;

%% Distribution of pre and post synaptic path lenghts by group

AlxPostSynapticPathLength = [];
AlxPreSynapticPathLength = [];
TransPostSynapticPathLength = []; 
TransPreSynapticPathLength = [];
DbxPostSynapticPathLength = [];
BarhlPostSynapticPathLength = [];


for i = 1:numel(cellIDs)
    if ismember(cellIDs{i}, cellIDsAlx) ==1
        AlxPostSynapticPathLength = [AlxPostSynapticPathLength; NormalizedLength(allLengthToPostNode{i})];
        AlxPreSynapticPathLength = [AlxPreSynapticPathLength; NormalizedLength(allLengthToPreNode{i})];
    elseif ismember(cellIDs{i}, cellIDsTrans) ==1
        TransPostSynapticPathLength = [TransPostSynapticPathLength; NormalizedLength(allLengthToPostNode{i})];
        TransPreSynapticPathLength = [TransPreSynapticPathLength ; NormalizedLength(allLengthToPreNode{i})];
    elseif ismember(cellIDs{i}, cellIDsDbx) ==1
        DbxPostSynapticPathLength = [DbxPostSynapticPathLength;NormalizedLength(allLengthToPostNode{i})];
    else
        BarhlPostSynapticPathLength = [BarhlPostSynapticPathLength;NormalizedLength(allLengthToPostNode{i})];
    end
end


subplot(4,1,1);
histogram(AlxPostSynapticPathLength,'FaceColor',[0.8,0,0],'BinWidth',0.1,'Normalization','probability'); % dimensions in microns
hold on;
line([mean(AlxPostSynapticPathLength),mean(AlxPostSynapticPathLength)], [0,0.35], 'Color',[0.9,0,0], 'LineWidth', 4);
histogram(AlxPreSynapticPathLength, 'FaceColor',[0,0.9,0],'BinWidth',0.1,'Normalization','probability'); % dimensions in microns
line([mean(AlxPreSynapticPathLength),mean(AlxPreSynapticPathLength)], [0,0.35], 'Color',[0,0.8,0], 'LineWidth', 4);

xlabel('Synaptic pathlenght in \mum', 'FontName', 'Arial', 'FontSize', 40);
ylabel('Probability', 'FontName', 'Arial', 'FontSize', 40);
set(gca, 'XLim',[0, 1],'YLim',[0,0.35],'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
box off;
set(gcf,'color','w');
set(gca,'color',[calx,0.2]);
daspect([1,2,1]);
% axis square;


subplot(4,1,2);
histogram(TransPostSynapticPathLength,'FaceColor',[0.8,0,0],'BinWidth',0.1,'Normalization','probability'); % dimensions in microns
hold on;
line([mean(TransPostSynapticPathLength),mean(TransPostSynapticPathLength)], [0,0.35], 'Color',[0.9,0,0], 'LineWidth', 4);
histogram(TransPreSynapticPathLength,'FaceColor',[0,0.9,0],'BinWidth',0.1,'Normalization','probability'); % dimensions in microns
line([mean(TransPreSynapticPathLength),mean(TransPreSynapticPathLength)], [0,0.35], 'Color',[0,0.8,0], 'LineWidth', 4);
%xlabel('Synaptic pathlenght in \mum', 'FontName', 'Arial', 'FontSize', 40);
%ylabel('Probability', 'FontName', 'Arial', 'FontSize', 40);
set(gca, 'XLim',[0, 1],'YLim',[0,0.35],'YColor','k','FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
box off;
set(gcf,'color','w');
set(gca,'color',[ctrans,0.2]);
daspect([1,2,1]);

% axis square;


subplot(4,1,3);
histogram(DbxPostSynapticPathLength,'FaceColor',[0.8,0,0],'BinWidth',0.1,'Normalization','probability'); % dimensions in microns
hold on;
line([mean(DbxPostSynapticPathLength),mean(DbxPostSynapticPathLength)], [0,0.35], 'Color',[0.9,0,0], 'LineWidth', 4);
%xlabel('Synaptic pathlenght in \mum', 'FontName', 'Arial', 'FontSize', 40);
%ylabel('Probability', 'FontName', 'Arial', 'FontSize', 40);
set(gca, 'XLim',[0, 1],'YLim',[0,0.35],'YColor','k','FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
box off;
set(gcf,'color','w');
set(gca,'color',[cdbx,0.2]);
daspect([1,2,1]);

%axis square;


subplot(4,1,4);
histogram(BarhlPostSynapticPathLength,'FaceColor',[0.8,0,0],'BinWidth',0.1,'Normalization','probability'); % dimensions in microns
hold on;
line([mean(BarhlPostSynapticPathLength),mean(BarhlPostSynapticPathLength)], [0,0.35], 'Color',[0.9,0,0], 'LineWidth', 4);
%xlabel('Synaptic pathlenght in \mum', 'FontName', 'Arial', 'FontSize', 40);
%ylabel('Probability', 'FontName', 'Arial', 'FontSize', 40);
set(gca, 'XLim',[0,1],'YLim',[0,0.35],'YColor','k','FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
box off;
set(gcf,'color','w');
set(gca,'color',[cbarhl,0.2]);
daspect([1,2,1]);
%axis square;


%% Minimum Distance between trees

index =0;
k=0;
clear minDistance;
clear minDistance1000;
figure();
for i = 1:size(cellIDs,2)
    if size(allPreSynapse{i},1)>0
        for j = 1:size(cellIDs,2)
            sz = 0;
            [MinDistance{i,j}] = MinSynapticDistance(allPreSynapse{i},allPost{j});
            sz = sz + sum(size(MinDistance{i,j},1)*size(MinDistance{i,j},2));
            minDistance(i,j) = min(min(MinDistance{i,j}));
            PairwiseSynapticDistance(index+1:index+sz,i) = [reshape(MinDistance{i,j},[],1)];
            index = sz;
        end
        clear temp;
        k=k+1;
        subplot(3,3,k);
        [a,b] = find(PairwiseSynapticDistance(:,i)>0 & PairwiseSynapticDistance(:,i)<5000);
        hist(PairwiseSynapticDistance(a,i)./1000); % report in microns
        xlim([0 5]);
        title(cellIDs{i});
        %figtitle('Synaptic pairs <5 \mum apart');
        clear a;
        clear b;
    else
        continue
    end
end
%export_fig('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/ZfishFigures/PairwiseSynapseDistrubution', '-eps');

% plot the minimum synaptic distance between trees
figure();
minDistance(all(~minDistance,2),:) = [];
imagesc(minDistance/1000); % in microns

colormap(parula);
c1 = colorbar;
c1.Label.String = 'Distance(\mum)';
%title('Minimum synaptic distance between trees in \mum');
xlabel('Postsynaptic cell', 'FontName', 'Arial', 'FontSize', 40);
ylabel('Presynaptic cell', 'FontName', 'Arial', 'FontSize', 40);
set(gca,'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
set(gcf,'color','w');
axis image;

%export_fig('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/ZfishFigures/MinimumSynapticDistance','-eps');

% plot only those synapses pairs that are within 1000nm of each other
figure();
[r,c] = find(minDistance>0 & minDistance<1000); % find synapses that are atleast 1000nm close to each other
minDistance1000 = zeros(size(cellIDs,2));
for i = 1:size(r,1)
    minDistance1000(r(i),c(i)) = minDistance(r(i),c(i));
end

minDistance1000(all(~minDistance1000,2),:) = [];

imagesc(minDistance1000/1000); % in microns
colormap(hot);
c2 = colorbar;
c2.Label.String = 'Distance(\mum)'
xlabel('Postsynaptic cell', 'FontName', 'Arial', 'FontSize', 40);
ylabel('Presynaptic cell', 'FontName', 'Arial', 'FontSize', 40);
set(gca, 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
set(gcf,'color','w');
axis image;
%export_fig('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/ZfishFigures/MinimumSynapticDistance1000', '-eps');

%% Euclidean Distances between soma and synaptic sites



AlxPreEucDist = [];
AlxPostEucDist = [];
TransPreEucDist = [];
TransPostEucDist = [];
DbxPostEucDist = [];
BarhlPostEucDist = [];


for i= 1:numel(cellIDs)
    i
    if ismember(cellIDs{i}, cellIDsAlx) ==1
        AlxPreEucDist =  [AlxPreEucDist; diag(pdist2(repmat(CellSoma(i,:),size(allPost{i},1),1),allPreSynapse{i}))];
        AlxPostEucDist = [AlxPostEucDist; diag(pdist2(repmat(CellSoma(i,:),size(allPost{i},1),1),allPost{i}))];
    elseif ismember(cellIDs{i}, cellIDsTrans) ==1
        TransPreEucDist =  [TransPreEucDist; diag(pdist2(repmat(CellSoma(i,:),size(allPost{i},1),1),allPreSynapse{i}))];
        TransPostEucDist = [TransPostEucDist; diag(pdist2(repmat(CellSoma(i,:),size(allPost{i},1),1),allPost{i}))];
    elseif ismember(cellIDs{i}, cellIDsDbx) ==1
        DbxPostEucDist = [DbxPostEucDist; diag(pdist2(repmat(CellSoma(i,:),size(allPost{i},1),1),allPost{i}))];
    else
        BarhlPostEucDist = [BarhlPostEucDist; diag(pdist2(repmat(CellSoma(i,:),size(allPost{i},1),1),allPost{i}))];
    end
end

figure();
subplot(2,2,1)
histogram(AlxPostEucDist/1000,'FaceColor',[0.8,0,0],'BinWidth',2,'Normalization','probability'); % dimensions in microns
hold on;
histogram(AlxPreEucDist/1000, 'FaceColor',[0,0.9,0],'BinWidth',2,'Normalization','probability'); % dimensions in microns
set(gca, 'XLim',[0, 120], 'YLim',[0,0.15],'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
box off;
set(gcf,'color','w');
set(gca,'color',[calx,0.2]);
axis square;

subplot(2,2,2)
histogram(TransPostEucDist/1000,'FaceColor',[0.8,0,0],'BinWidth',2,'Normalization','probability'); % dimensions in microns
hold on;
histogram(TransPreEucDist/1000,'FaceColor',[0,0.9,0],'BinWidth',2,'Normalization','probability'); % dimensions in microns
set(gca, 'XLim',[0, 120],'YLim',[0,0.15],'XTickLabel',[],'YColor','k','FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
box off;
set(gcf,'color','w');
set(gca,'color',[ctrans,0.2]);
axis square;

subplot(2,2,3)
histogram(DbxPostEucDist/1000,'FaceColor',[0.8,0,0],'BinWidth',2,'Normalization','probability'); % dimensions in microns
set(gca, 'XLim',[0, 120],'YLim',[0,0.15],'XTickLabel',[],'YColor','k','FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
box off;
set(gcf,'color','w');
set(gca,'color',[cdbx,0.2]);
axis square;

subplot(2,2,4)
histogram(BarhlPostEucDist/1000,'FaceColor',[0.8,0,0],'BinWidth',2,'Normalization','probability'); % dimensions in microns
set(gca, 'XLim',[0, 120],'YLim',[0,0.15],'XTickLabel',[],'YColor','k','FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
box off;
set(gcf,'color','w');
set(gca,'color',[cbarhl,0.2]);
axis square;
xlabel('Euc. distance (\mum)', 'FontName', 'Arial', 'FontSize', 40);
ylabel('Probability', 'FontName', 'Arial', 'FontSize', 40);





