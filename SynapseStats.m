% Synapse Stats

calx = [1 0.5 0.3];
cdbx = [1 0.3 1];
cbarhl = [0.3 0.5 1];
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
hold on;



[ax,h1,h2] = plotyy([ones(length(AlxPost),1)' , 2*ones(length(DbxPost),1)',3*ones(length(BarhlPost),1)'] , [AlxPost,DbxPost,BarhlPost], 1.2*ones(length(AlxPre),1)', AlxPre) ; % 'o', );
h1.Marker = 'o';
h1.MarkerFaceColor =  [0.7,0,0];
h1.MarkerEdgeColor = 'none';
h1.MarkerSize =  25 ;
h1.LineStyle = 'none';
h2.Marker = 'o';
h2.MarkerFaceColor =  [0,0.7,0];
h2.MarkerEdgeColor = 'none';
h2.MarkerSize =  25;
h2.LineStyle = 'none';


h3 = plot(1.2,MeanPre,'o', 'MarkerFaceColor', calx, 'MarkerEdgeColor','none', 'MarkerSize', 35 );
h4 = plot([1.2;1.2], [MeanPre-StdPre;MeanPre+StdPre], 'Color','k','LineWidth',2);
%uistack(h4,'top');

plot(1,MeanPost(1), 'o', 'MarkerFaceColor', calx, 'MarkerEdgeColor','none', 'MarkerSize', 35 );
plot(2,MeanPost(2), 'o', 'MarkerFaceColor', cdbx, 'MarkerEdgeColor','none', 'MarkerSize', 35 );
plot(3,MeanPost(3), 'o', 'MarkerFaceColor', cbarhl, 'MarkerEdgeColor','none', 'MarkerSize', 35 );
plot([1:3;1:3], [MeanPost-StdPost;MeanPost+StdPost], 'Color','k','LineWidth',2);

set(ax(1),'XTick', [1:3],'XTickLabel', {'group1'; 'group2'; 'group3'},'XLim', [0.5 3.5], 'FontName', 'Arial', 'FontSize', 40, 'Ycolor', [0.7,0,0], 'LineWidth',2);
set(ax(2),'XTick', [1:3],'XTickLabel', {'group1'; 'group2'; 'group3'}, 'XLim', [0.5 3.5],'FontName', 'Arial', 'FontSize', 40, 'Ycolor' , [0,0.7,0],'LineWidth',2);
set(h1,'color','w');
set(h2,'color','w');
axis (ax(1),'square');
axis(ax(2),'square');
ylabel(ax(1),'Number of postsynaptic sites', 'FontName', 'Arial', 'FontSize', 40);
ylabel(ax(2), 'Number of presynaptic sites', 'FontName', 'Arial', 'FontSize', 40);
set(gcf,'color','w');
box off;

%% Distribution of inter-synaptic distance

clear PostSynapticDistance;
clear tempPost;
PostSynapticDistance = [];
figure;

for ii = 1:size(cellIDs,2)
    [y,I] = sort(allLengthToPostNode{ii});
    for i = 2:length(I)
        tempPost(ii,i-1) = abs(allLengthToPostNode{ii}(I(i-1)) - allLengthToPostNode{ii}(I(i)));
    end
    PostSynapticDistance = [PostSynapticDistance;tempPost(ii,:)'];
end
histogram(find(PostSynapticDistance>0)/1000, 'FaceColor',[0.5,0.5,0.5]);
%title('Inter-postsynaptic distance', 'FontName', 'Arial', 'FontSize', 20);
xlabel('Inter postsynaptic distance in \mum', 'FontName', 'Arial', 'FontSize', 40);
ylabel('Count', 'FontName', 'Arial', 'FontSize', 40);
set(gca, 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
box off;set(gcf,'color','w');
axis square;



clear PreSynapticDistance;
clear tempPre;
PreSynapticDistance = [];
tempPre =[];
figure();
index =1;
for ii = 1:size(cellIDs,2)
    if ~isempty(allLengthToPreNode{ii})==1 & size(allLengthToPreNode{ii},1)>1
        [y,I] = sort(allLengthToPreNode{ii});
        for i = 2:length(I)
            tempPre(ii,i-1) = abs(allLengthToPreNode{ii}(I(i-1)) - allLengthToPreNode{ii}(I(i)));
        end
        PreSynapticDistance = [PreSynapticDistance;tempPre(index,:)'];
        index = index+1;
    else
        index = index+1;
        continue;
    end
end

histogram(find(PreSynapticDistance>0)/1000, 'FaceColor',[0.5,0.5,0.5]);
%title('Inter-presynaptic distance');
xlabel('Inter presynaptic distance in \mum', 'FontName', 'Arial', 'FontSize', 40);
ylabel('Count', 'FontName', 'Arial', 'FontSize', 40);
set(gca, 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
box off;
set(gcf,'color','w');
axis square
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
histogram(allPostSynapticLength/1000,'FaceColor',[0.5,0.5,0.5]); % dimensions in microns
xlabel('Postsynaptic pathlenght in \mum', 'FontName', 'Arial', 'FontSize', 40);
ylabel('Count', 'FontName', 'Arial', 'FontSize', 40);
set(gca, 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
box off;
set(gcf,'color','w');
axis square

figure();
histogram(allPreSynapticLength/1000,'FaceColor',[0.5,0.5,0.5]); % dimensions in microns
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

