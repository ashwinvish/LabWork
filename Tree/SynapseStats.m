% Synapse Stats


AlxPre = [];
AlxPost = [];
DbxPost = [];
BarhlPost = [];

for kk = 1:numel(cellIDs)
    % subplot(3,8,kk);
    if ismember(cellIDs{kk},cellIDsAlx)==1
        AlxPost  = [AlxPost, length(allPost{kk})];
        AlxPre = [AlxPre, length(allPreSynapse{kk})];
    elseif ismember(cellIDs{kk},cellIDsDbx)==1
        DbxPost = [DbxPost, length(allPost{kk})];
    else
        BarhlPost = [BarhlPost, length(allPost{kk})];
    end
end

% bar plot of number of Postsynapses
figure();
MeanPost = [mean(AlxPost), mean(DbxPost), mean(BarhlPost)];
StdPost = [std(AlxPost), std(DbxPost), std(BarhlPost)];
MeanPre = mean(AlxPre);
StdPre = std(AlxPre);


plot([ones(length(AlxPost),1)' , 2*ones(length(DbxPost),1)',3*ones(length(BarhlPost),1)'],[AlxPost,DbxPost,BarhlPost], 'Marker', 'o', 'MarkerFaceColor', [0.7,0,0],'MarkerSize',25,'LineStyle','none', 'MarkerEdgeColor','k');
hold on;
plot(1.2*ones(length(AlxPre),1)',AlxPre, 'Marker', 'o', 'MarkerFaceColor', [0,0.7,0],'MarkerSize',25,'LineStyle','none',  'MarkerEdgeColor','k');
plot(1.2,MeanPre,'o', 'MarkerFaceColor', calx, 'MarkerEdgeColor','k', 'MarkerSize', 35 );
plot([1.2;1.2], [MeanPre-StdPre;MeanPre+StdPre], 'Color','k','LineWidth',2);
plot(1,MeanPost(1), 'o', 'MarkerFaceColor', calx, 'MarkerEdgeColor','k', 'MarkerSize', 35 );
plot(2,MeanPost(2), 'o', 'MarkerFaceColor', cdbx, 'MarkerEdgeColor','k', 'MarkerSize', 35 );
plot(3,MeanPost(3), 'o', 'MarkerFaceColor', cbarhl, 'MarkerEdgeColor','k', 'MarkerSize', 35 );
plot([1:3;1:3], [MeanPost-StdPost;MeanPost+StdPost], 'Color','k','LineWidth',2);

set(gca, 'XTick', [1:3],'XTickLabel', {'group1'; 'group2'; 'group3'},'XLim', [0.5 3.5], 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
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

histogram(InterPostSynapticDistance/1000, 'FaceColor',[0.5,0.5,0.5],'BinWidth',0.5);
xlabel('Inter postsynaptic distance in \mum', 'FontName', 'Arial', 'FontSize', 40);
ylabel('Count', 'FontName', 'Arial', 'FontSize', 40);
set(gca,'XTick',[0,5,10,15],'XLim',[-0.5,15], 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
box off;set(gcf,'color','w');
axis square;

clear InterPreSynapticDistance;
InterPreSynapticDistance = [];
tempPre =[];
figure();

for ii = 1:size(cellIDs,2)
    if ~isempty(allPreDiff{ii})
        InterPreSynapticDistance = [InterPreSynapticDistance; nonzeros(allPreDiff{ii})];
    end
end

histogram(InterPreSynapticDistance/1000, 'FaceColor',[0.5,0.5,0.5],'BinWidth',0.5);
xlabel('Inter presynaptic distance in \mum', 'FontName', 'Arial', 'FontSize', 40);
ylabel('Count', 'FontName', 'Arial', 'FontSize', 40);
set(gca, 'XTick',[0,5,10,15,20,25],'XLim',[-0.5,25],'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
box off;
set(gcf,'color','w');
axis square

figure();
histogram(InterPostSynapticDistance/1000, 'FaceColor',[0.9,0.0,0],'BinWidth',0.5, 'Normalization','probability');
hold on;
histogram(InterPreSynapticDistance/1000, 'FaceColor',[0,0.8,0],'BinWidth',0.5, 'Normalization','probability');
axis square;
xlabel('Inter synaptic distance in \mum', 'FontName', 'Arial', 'FontSize', 40);
ylabel('Probability', 'FontName', 'Arial', 'FontSize', 40);
set(gca, 'XTick',[0,5,10,15,20,25],'XLim',[-0.5,25],'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
box off;
set(gcf,'color','w');

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
histogram(allPostSynapticLength/1000,'FaceColor',[0.5,0.5,0.5],'BinWidth',10); % dimensions in microns
xlabel('Postsynaptic pathlenght in \mum', 'FontName', 'Arial', 'FontSize', 40);
ylabel('Count', 'FontName', 'Arial', 'FontSize', 40);
set(gca, 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
box off;
set(gcf,'color','w');
axis square

figure();
histogram(allPreSynapticLength/1000,'FaceColor',[0.5,0.5,0.5],'BinWidth',10); % dimensions in microns
%title('Distribution of Presynaptic pathlength');
xlabel('Presynaptic pathlenght in \mum', 'FontName', 'Arial', 'FontSize', 40);
ylabel('Count', 'FontName', 'Arial', 'FontSize', 40);
set(gca, 'XLim',[0, 250],'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
box off;
set(gcf,'color','w');
axis square

%% Minimum Distance between trees

index =0;
k=0;
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
        figtitle('Synaptic pairs <5 \mum apart');
        clear a;
        clear b;
    else
        continue
    end
end
%export_fig('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/ZfishFigures/PairwiseSynapseDistrubution', '-eps');

% plot the minimum synaptic distance between trees
figure();
imagesc(minDistance/1000); % in microns
colorbar;
%title('Minimum synaptic distance between trees in \mum');
xlabel('Postsynaptic cell', 'FontName', 'Arial', 'FontSize', 40);
ylabel('Presynaptic cell', 'FontName', 'Arial', 'FontSize', 40);
set(gca, 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
box off;
set(gcf,'color','w');
axis square;

%export_fig('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/ZfishFigures/MinimumSynapticDistance','-eps');

% plot only those synapses pairs that are within 1000nm of each other
figure();
[r,c] = find(minDistance>0 & minDistance<1000); % find synapses that are atleast 1000nm close to each other
minDistance1000 = zeros(size(cellIDs,2));
for i = 1:size(r,1)
    minDistance1000(r(i),c(i)) = minDistance(r(i),c(i));
end
imagesc(minDistance1000/1000); % in microns
colorbar;

xlabel('Postsynaptic cell', 'FontName', 'Arial', 'FontSize', 40);
ylabel('Presynaptic cell', 'FontName', 'Arial', 'FontSize', 40);
set(gca, 'FontName', 'Arial', 'FontSize', 40), 'LineWidth',2;
box off;
set(gcf,'color','w');
axis square;
%export_fig('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/ZfishFigures/MinimumSynapticDistance1000', '-eps');

