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
        BarhlPost = [DbxPost, length(allPost{kk})];
    end
end

% bar plot of number of Postsynapses
figure();
MeanPost = [mean(AlxPost), mean(DbxPost), mean(BarhlPost)];
StdPost = [std(AlxPost), std(DbxPost), std(BarhlPost)];
MeanPre = mean(AlxPre);
StdPre = std(AlxPre);

[ax,h1,h2] = plotyy([ones(length(AlxPost),1)' , 2*ones(length(DbxPost),1)',3*ones(length(BarhlPost),1)'] , [AlxPost,DbxPost,BarhlPost], 1.2*ones(length(AlxPost),1)', AlxPre) ; % 'o', );
h1.Marker = 'o';
h1.MarkerFaceColor =  [0.7,0.7,0.7];
h1.MarkerEdgeColor = 'none';
h1.MarkerSize =  25 ;
h1.LineStyle = 'none';
h2.Marker = 'o';
h2.MarkerFaceColor =  [0.9,0.9,0.9];
h2.MarkerEdgeColor = 'none';
h2.MarkerSize =  25;
h2.LineStyle = 'none';
hold on;
plot(1,MeanPost(1), 'o', 'MarkerFaceColor', calx, 'MarkerEdgeColor','none', 'MarkerSize', 25 );
plot(2,MeanPost(2), 'o', 'MarkerFaceColor', cdbx, 'MarkerEdgeColor','none', 'MarkerSize', 25 );
plot(3,MeanPost(3), 'o', 'MarkerFaceColor', cbarhl, 'MarkerEdgeColor','none', 'MarkerSize', 25 );
plot([1:3;1:3], [MeanPost-StdPost;MeanPost+StdPost], 'Color','k','LineWidth',2);

plot(1.2,MeanPre,'o', 'MarkerFaceColor', calx, 'MarkerEdgeColor','none', 'MarkerSize', 25 );
h3 = plot([1.2;1.2], [MeanPre-StdPre;MeanPre+StdPre], 'Color','k','LineWidth',2);
uistack(h3,'top');

set(ax(1),'XTick', [1:3],'XTickLabel', {'group1'; 'group2'; 'group3'},'XLim', [0.5 3.5], 'FontName', 'Arial', 'FontSize', 20, 'Ycolor','k');
set(ax(2),'XTick', [1:3],'XTickLabel', {'group1'; 'group2'; 'group3'}, 'XLim', [0.5 3.5],'FontName', 'Arial', 'FontSize', 20, 'Ycolor' , 'k');


set(h1,'color','w');
set(h2,'color','w');

axis (ax(1),'square');
axis(ax(2),'square');

ylabel(ax(1),'Number of postsynaptic sites', 'FontName', 'Arial', 'FontSize', 20);
ylabel(ax(2), 'Number of presynaptic sites', 'FontName', 'Arial', 'FontSize', 20);
box off;

%% Distribution of inter-synaptic distance

clear PostSynapticDistance;
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
xlabel('Distance in \mum', 'FontName', 'Arial', 'FontSize', 20);
ylabel('Count', 'FontName', 'Arial', 'FontSize', 20);
set(gca, 'FontName', 'Arial', 'FontSize', 20);
box off;set(gcf,'color','w');
axis square;



clear PreSynapticDistance;
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
xlabel('Distance in \mum', 'FontName', 'Arial', 'FontSize', 20);
ylabel('Count', 'FontName', 'Arial', 'FontSize', 20);
set(gca, 'FontName', 'Arial', 'FontSize', 20);
box off;
set(gcf,'color','w');
axis square

%%

