% Synapse Stats


AlxPre = [];
AlxPost = [];
TransPre = [];
TransPost = [];
DbxPost = [];
BarhlPost = [];

for kk = 1:numel(cellIDs)
    % subplot(3,8,kk);
    if ismember(cellIDs{kk},cellIDsAlx)==1
        AlxPost  = [AlxPost, length(allPost{kk})];
        AlxPre = [AlxPre, length(allPreSynapse{kk})];
    elseif ismember(cellIDs{kk},cellIDsTrans)==1
        TransPost  = [TransPost, length(allPost{kk})];
        TransPre = [TransPre, length(allPreSynapse{kk})];
    elseif ismember(cellIDs{kk},cellIDsDbx)==1
        DbxPost = [DbxPost, length(allPost{kk})];
    else
        BarhlPost = [BarhlPost, length(allPost{kk})];
    end
end

% bar plot of number of Postsynapses
figure();
MeanPost = [mean(AlxPost), mean(TransPost), mean(DbxPost), mean(BarhlPost)];
StdPost = [std(AlxPost), std(TransPost), std(DbxPost), std(BarhlPost)];
MeanPre = [mean(AlxPre), mean(TransPre)];
StdPre = [std(AlxPre), std(TransPre)];


plot([ones(length(AlxPost),1)',2*ones(length(TransPost),1)', 3*ones(length(DbxPost),1)',4*ones(length(BarhlPost),1)'],[AlxPost,TransPost,DbxPost,BarhlPost], 'Marker', 'o', 'MarkerFaceColor', [0.7,0,0],'MarkerSize',25,'LineStyle','none', 'MarkerEdgeColor','k');
hold on;
plot(1.2*ones(length(AlxPre),1)',AlxPre, 'Marker', 'o', 'MarkerFaceColor', [0,0.7,0],'MarkerSize',25,'LineStyle','none',  'MarkerEdgeColor','k');
plot(2.2*ones(length(TransPre),1)',TransPre, 'Marker', 'o', 'MarkerFaceColor', [0,0.7,0],'MarkerSize',25,'LineStyle','none',  'MarkerEdgeColor','k' );

plot(1.2,MeanPre(1),'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor','k', 'MarkerSize', 25 );
plot(2.2,MeanPre(2),'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor','k', 'MarkerSize', 25 );

plot([1.2,2.2;1.2,2.2], [MeanPre-StdPre;MeanPre+StdPre], 'Color','k','LineWidth',2);
plot(1,MeanPost(1), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor','k', 'MarkerSize', 25 );
plot(2,MeanPost(2), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor','k', 'MarkerSize', 25 );
plot(3,MeanPost(3), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor','k', 'MarkerSize', 25 );
plot(4,MeanPost(4), 'o', 'MarkerFaceColor',' k', 'MarkerEdgeColor','k', 'MarkerSize', 25 );
plot([1:4;1:4], [MeanPost-StdPost;MeanPost+StdPost], 'Color','k','LineWidth',2);

set(gca, 'XTick', [1:4],'XTickLabel', {'group1'; 'group2'; 'group3'; 'group4'},'XLim', [0.5 4.5], 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
ylabel('Number of sites', 'FontName', 'Arial', 'FontSize', 40);
set(gcf,'color','w');
legend('Postsynapse','Presynapse');
legend('boxoff');
box off;
axis square;

%% Distribution of inter-synaptic distance

clear InterPostSynapticDistance;
InterPostSynapticDistance = [];
figure;

for ii = 1:numel(cellIDs)
    InterPostSynapticDistance = [InterPostSynapticDistance;nonzeros(allPostDiff{ii})];
end

histogram(InterPostSynapticDistance/1000, 'FaceColor',[0.9,0,0],'BinWidth',500/1000, 'Normalization','probability');
hold on;
xlabel('Inter postsynaptic distance (\mum)', 'FontName', 'Arial', 'FontSize', 40);
ylabel('Proabability', 'FontName', 'Arial', 'FontSize', 40);
set(gca,'XTick',[0,5,10,15],'XLim',[-0.5,15], 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
axis square;
% lognormal distribution
pd1 = fitdist(InterPostSynapticDistance/1000,'lognormal');
x_val =  min(InterPostSynapticDistance/1000):10/1000:max(InterPostSynapticDistance/1000);
y_val = pdf(pd1,x_val);
plot(x_val,y_val,'k','LineWidth',2);
box off;
set(gcf,'color','w');

clear x_val;
clear y_val;

clear InterPreSynapticDistance;
InterPreSynapticDistance = [];
tempPre =[];
figure();

for ii = 1:size(cellIDs,2)
    if ~isempty(allPreDiff{ii})
        InterPreSynapticDistance = [InterPreSynapticDistance; nonzeros(allPreDiff{ii})];
    end
end

histogram(InterPreSynapticDistance/1000, 'FaceColor',[0,0.8,0],'BinWidth',0.5, 'Normalization','Probability');
hold on;
xlabel('Inter presynaptic distance (\mum)', 'FontName', 'Arial', 'FontSize', 40);
ylabel('Probability', 'FontName', 'Arial', 'FontSize', 40);
set(gca, 'XTick',[0,5,10,15,20,25],'XLim',[-0.5,25],'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);


%lognormal  distribution
pd2 = fitdist(InterPreSynapticDistance/1000,'lognormal');
x_val =  min(InterPreSynapticDistance/1000):10/1000:max(InterPreSynapticDistance/1000);
y_val = pdf(pd2,x_val);
plot(x_val,y_val,'k','LineWidth',2);
box off;
set(gcf,'color','w');


figure();
histogram(InterPostSynapticDistance/1000, 'FaceColor',[0.9,0.0,0],'BinWidth',0.5, 'Normalization','probability');
hold on;
histogram(InterPreSynapticDistance/1000, 'FaceColor',[0,0.8,0],'BinWidth',0.5, 'Normalization','probability');
axis square;
xlabel('Inter synaptic distance in \mum', 'FontName', 'Arial', 'FontSize', 40);
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


figure() 
histogram(AlxInterPost/1000, 'FaceColor',[0.9,0,0], 'BinWidth',0.5, 'Normalization','probability');
hold on;
histogram(AlxInterPre/1000, 'FaceColor',[0,0.8,0],'BinWidth',0.5,'Normalization','probability');
xlabel('Synaptic pathlenght (\mum)', 'FontName', 'Arial', 'FontSize', 40);
ylabel('Probability', 'FontName', 'Arial', 'FontSize', 40);
set(gca, 'XLim',[0, 25], 'YLim',[0,0.4],'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
box off;
set(gcf,'color','w');
set(gca,'color',[calx,0.2]);
axis square;


figure();
histogram(TransInterPost/1000, 'FaceColor',[0.9,0,0], 'BinWidth',0.5, 'Normalization','probability');
hold on;
histogram(TransInterPre/1000,'FaceColor',[0,0.8,0],'BinWidth',0.5,'Normalization','probability');
ylabel('Probability', 'FontName', 'Arial', 'FontSize', 40);
set(gca, 'XLim',[0, 25],'YLim',[0,0.4],'XTickLabel',[],'YColor','none','FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
box off;
set(gcf,'color','w');
set(gca,'color',[ctrans,0.2]);
axis square;

figure();
histogram(DbxInterPost/1000,  'FaceColor',[0.9,0,0], 'BinWidth',0.5, 'Normalization','probability');
ylabel('Probability', 'FontName', 'Arial', 'FontSize', 40);
set(gca, 'XLim',[0, 25],'YLim',[0,0.4],'XTickLabel',[],'YColor','none','FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
box off;
set(gcf,'color','w');
set(gca,'color',[cdbx,0.2]);
axis square;

figure();
histogram(BarhlInterPost/1000, 'FaceColor',[0.9,0,0], 'BinWidth',0.5, 'Normalization','probability');
ylabel('Probability', 'FontName', 'Arial', 'FontSize', 40);
set(gca, 'XLim',[0, 25],'YLim',[0,0.4],'XTickLabel',[],'YColor','none','FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
box off;
set(gcf,'color','w');
set(gca,'color',[cbarhl,0.2]);
axis square;







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
histogram(allPostSynapticLength/1000,'FaceColor',[0.8,0,0],'BinWidth',50,'Normalization','probability'); % dimensions in microns
%xlabel('Postsynaptic pathlenght in \mum', 'FontName', 'Arial', 'FontSize', 40);
ylabel('Count', 'FontName', 'Arial', 'FontSize', 40);
set(gca, 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
box off;
set(gcf,'color','w');
axis square

hold on;


histogram(allPreSynapticLength/1000,'FaceColor',[0,0.9,0],'BinWidth',50,'Normalization','probability'); % dimensions in microns
%title('Distribution of Presynaptic pathlength');
xlabel('Synaptic pathlenght in \mum', 'FontName', 'Arial', 'FontSize', 40);
ylabel('Probability', 'FontName', 'Arial', 'FontSize', 40);
legend({'PostSynapse','PreSynapse'}, 'box','off');
set(gca,'XLim',[0,1250],'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
box off;
set(gcf,'color','w');
axis square

%% Distribution of pre and post synaptic path lenghts by group

allAlxPostSynapticPathLength = [];
allAlxPreSynapticPathLength = [];
allTransPostSynapticPathLength = []; 
allTransPreSynapticPathLength = [];
allDbxPostSynapticPathLength = [];
allBarhlPostSynapticPathLength = [];


for i = 1:numel(cellIDs)
    if ismember(cellIDs{i}, cellIDsAlx) ==1
        allAlxPostSynapticPathLength = [allAlxPostSynapticPathLength; allLengthToPostNode{i}];
        allAlxPreSynapticPathLength = [allAlxPreSynapticPathLength; allLengthToPreNode{i}];
    elseif ismember(cellIDs{i}, cellIDsTrans) ==1
        allTransPostSynapticPathLength = [  allTransPostSynapticPathLength; allLengthToPostNode{i}];
        allTransPreSynapticPathLength = [allTransPreSynapticPathLength ; allLengthToPreNode{i}];
    elseif ismember(cellIDs{i}, cellIDsDbx) ==1
        allDbxPostSynapticPathLength = [allDbxPostSynapticPathLength;allLengthToPostNode{i}];
    else
        allBarhlPostSynapticPathLength = [ allBarhlPostSynapticPathLength;allLengthToPostNode{i}];
    end
end


figure();
histogram(allAlxPostSynapticPathLength/1000, 'DisplayStyle','stairs','EdgeColor',[0.8,0,0],'BinWidth',20,'Normalization','probability','LineWidth',2); % dimensions in microns
hold on;
histogram(allAlxPreSynapticPathLength/1000, 'DisplayStyle','stairs','EdgeColor',[0,0.9,0],'BinWidth',20,'Normalization','probability','LineWidth',2); % dimensions in microns
xlabel('Synaptic pathlenght in \mum', 'FontName', 'Arial', 'FontSize', 40);
ylabel('Probability', 'FontName', 'Arial', 'FontSize', 40);
set(gca, 'XLim',[0, 1200], 'YLim',[0,0.15],'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
box off;
set(gcf,'color','w');
set(gca,'color',[calx,0.2]);
axis square;

figure();
histogram(allTransPostSynapticPathLength/1000,'DisplayStyle','stairs','EdgeColor',[0.8,0,0],'BinWidth',20,'Normalization','probability','LineWidth',2); % dimensions in microns
hold on;
histogram(allTransPreSynapticPathLength/1000,'DisplayStyle','stairs','EdgeColor',[0,0.9,0],'BinWidth',20,'Normalization','probability','LineWidth',2); % dimensions in microns
%xlabel('Synaptic pathlenght in \mum', 'FontName', 'Arial', 'FontSize', 40);
%ylabel('Probability', 'FontName', 'Arial', 'FontSize', 40);
set(gca, 'XLim',[0, 1200],'YLim',[0,0.15],'XTickLabel',[],'YColor','none','FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
box off;
set(gcf,'color','w');
set(gca,'color',[ctrans,0.2]);
axis square;

figure();
histogram(allDbxPostSynapticPathLength/1000,'DisplayStyle','stairs','EdgeColor',[0.8,0,0],'BinWidth',20,'Normalization','probability','LineWidth',2); % dimensions in microns
%xlabel('Synaptic pathlenght in \mum', 'FontName', 'Arial', 'FontSize', 40);
%ylabel('Probability', 'FontName', 'Arial', 'FontSize', 40);
set(gca, 'XLim',[0, 1200],'YLim',[0,0.15],'XTickLabel',[],'YColor','none','FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
box off;
set(gcf,'color','w');
set(gca,'color',[cdbx,0.2]);
axis square;

figure();
histogram(allBarhlPostSynapticPathLength/1000,'DisplayStyle','stairs','EdgeColor',[0.8,0,0],'BinWidth',20,'Normalization','probability','LineWidth',2); % dimensions in microns
%xlabel('Synaptic pathlenght in \mum', 'FontName', 'Arial', 'FontSize', 40);
%ylabel('Probability', 'FontName', 'Arial', 'FontSize', 40);
set(gca, 'XLim',[0, 1200],'YLim',[0,0.15],'XTickLabel',[],'YColor','none','FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
box off;
set(gcf,'color','w');
set(gca,'color',[cbarhl,0.2]);
axis square;


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
set(gca, 'FontName', 'Arial', 'FontSize', 40), 'LineWidth',2;
set(gcf,'color','w');
axis image;
%export_fig('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/ZfishFigures/MinimumSynapticDistance1000', '-eps');

%% Some stats


