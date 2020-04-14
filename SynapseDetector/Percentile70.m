clc
clear;

load AllCells.mat
load ConnMatrixPre.mat
LoadDataFrame;

[m,n,v] = find(ConnMatrixPre);           % m -postSynapse,row, n - column, v - number of synapses,value
G = digraph(n',m',v',size(AllCells,1)); % digraph(source, target, weight, numNodes)
G.Nodes.Names = num2str(AllCells);      % nodeNames

for i = 1:length(AllCells)
temp = SynapticPartners(AllCells(i),1,df);
inDegVect(i) = length(temp<1e5); % number of presynapses
faninVect(i) = length(unique(temp(temp<1e5))); % number of presynaptic partners
clear temp;
end

% potential synapses
allAxonsSites = [];
for i = 1:length(AllCells)
allAxonsSites = [allAxonsSites; df.presyn_x(df.postsyn_segid == AllCells(i)),...
                  df.presyn_y(df.postsyn_segid == AllCells(i)),...
                  df.presyn_z(df.postsyn_segid == AllCells(i))];        
end
allAxonsSites = TransformPoints(allAxonsSites,0);

for i = 1:length(AllCells)
    postSynapticSites = [df.postsyn_x(df.postsyn_segid == AllCells(i)),...
                  df.postsyn_y(df.postsyn_segid == AllCells(i)),...
                  df.postsyn_z(df.postsyn_segid == AllCells(i))];
      postSynapticSites = TransformPoints(postSynapticSites,0);        
   for jj = 1:size(postSynapticSites,1)        
             distmat =  pdist2(postSynapticSites(jj,:),allAxonsSites);
             potSites_2{i}(jj) = length(find(distmat<2));
             potSites_5{i}(jj) = length(find(distmat<5));
             potSites_10{i}(jj) = length(find(distmat<10));
             clear distmat;
   end
   clear postSynapticSites
end
%NodeForClustering = G.Nodes(find(G.indegree>prctile(G.indegree,70) & G.outdegree>prctile(G.outdegree,70)),1);
figure;

subplot(441); loglog(inDegVect, cellfun(@sum,potSites_2),'.');axis square; title('2um');
box off;

subplot(442); loglog(inDegVect, cellfun(@sum,potSites_5),'.');axis square; title('5um');
box off;

subplot(443); loglog(inDegVect, cellfun(@sum,potSites_10),'.');axis square; title('10um');
box off;
xlabel('in degree')
ylabel('potential synapses')


clear m
clear n
clear v
%%
index = find(ismember(AllCells,str2num(NodeForClustering.Names) ));
connMat_syn = ConnMatrixPre(index,index);
[m,n,v] = find(connMat_syn);                    % m -postSynapse,row, n - column, v - number of synapses,value
AllCellNames = AllCells(index);
% 
G1 = digraph(n',m',v',size(AllCellNames,1)); % digraph(source, target, weight, numNodes)
G1.Nodes.Names = num2str(AllCellNames);      % nodeNames
nn = numnodes(G1);
[s,t] = findedge(G1);
A = sparse(s,t,G1.Edges.Weight,nn,nn);


 gamma = 0.1:0.1:2;
 meanMod = [];
 meanClust = [];
% Modu = [];
%meanRandMod = [];
     jj = 1
 for kk = 1:10
     randMat = randomize_graph_partial_und(A,randi([0,1],size(A)),1000); %
%      for jj = 1:length(gamma)
        Q0 = -1; Q1 = 0;                            % initialize modularity values
        i =1 ;
        while abs(Q1-Q0)>1e-5                        % while modularity increases
            Q0 = Q1;                                  % perform community detection
            [M, Q1] = community_louvain(A,1,1:size(connMat_syn,1));
            %[M2,Q2] = community_louvain(randMat,gamma(jj),1:size(randMat,1));
            i = i+1;
        end
         meanMod(kk,jj) = Q1;
         cluster(:,kk) = M;
         meanClust(kk,jj) = length(unique(M));
        % meanRandMod(kk,jj) = Q2;
     end
% end
% figure;
% shadedErrorBar(gamma,mean(meanMod-meanRandMod),std(meanMod-meanRandMod)); axis square;


% association matrix
D = agreement(cluster); 
% random association matrix;
Tr = randi([min(unique(cluster)),max(unique(cluster))],[273,kk]);
Dr = agreement(Tr);
% thresholded agreement matrix;
T = consensus_und(D,max(Dr(:)),100);

%%

[Vsort,idx] = sort(T,'descend');
 
%subplot(2,2,1)
tempMat = connMat_syn(idx,idx);
tempMapThreshold = tempMat;
figure;
subplot(4,4,1)
histogram(tempMapThreshold,1:22,'EdgeColor','k','LineWidth',2,'FaceColor','w');
box off;
axis square;
set(gca,'YScale','log');
xlabel('Synapse count');
ylabel('Number of connections');


figure(2); 
tempMapThreshold(tempMapThreshold(:)>5) = 5;
cspy(tempMapThreshold,'Colormap',colorcet('R3','N',5),'Levels',5,'MarkerSize',10);
blocks = find(diff(round(Vsort)));
c = colorbar;
 c.Ticks = [1:5];
 c.TickLabels = [1:5];
hold on;

for i = 1:length(blocks)
    line([0,size(connMat_syn,2)]+0.5,[blocks(i),blocks(i)]+0.5,'color','k');
    line([blocks(i),blocks(i)]+0.5,[0,size(connMat_syn,2)]+0.5,'color','k');
end
box on;
%set(gca,'Xtick',[],'Ytick',[]);

%%

ABDr_CellIDs = [82140,82145,82143,77648,82146,77302,77301,82194,77705,77710,77305,77672,77300,77709,77661];
ABDc_CellIDs = [82213,77654,77646,82212,77295,81172,77657,82197,77682,77154, 77652,77658,77628,77292,77688,82195,77296];
ABDIr_CellIDs = [78553,77150,77886,78547,77631,77158,77665,78552,77668,77618,77634,78556];
ABDIc_CellIDs = [78574,79051,77148,77625,77144,77692,77641,79066,77640,77643];

Rov3 = [79961,81410,77449,83151,77456,77441,81611];
MiV2 = [78577,77694,77695,78566];
MiV1 = [80327,82221,82220,78940,82218,82217];
MiM1 = [77265];
RoM2 = [79244];
Mid3i = [76562];
RoL1 = [79395];
RSall = [Rov3,MiV2,MiV1,MiM1,RoM2,Mid3i,RoL1];


% 
% matOrder_commdet = [block1.cellIDs', block2.cellIDs',...
%     ABDr_CellIDs,ABDc_CellIDs,ABDIr_CellIDs,ABDIc_CellIDs];
clear matOrder_commdet_subMod;
matOrder_commdet_subMod = [AllCellNames(idx)',...
    ABDr_CellIDs,ABDc_CellIDs,ABDIr_CellIDs,ABDIc_CellIDs];
clear matIndex_commdet
[~,matIndex_commdet] = ismember(matOrder_commdet_subMod,AllCells);

gap = 20;
tempMat = zeros(length(matOrder_commdet_subMod)+gap);

tempMat(1:length(AllCellNames),1:length(AllCellNames)) = ConnMatrixPre(matIndex_commdet(1:length(AllCellNames)),matIndex_commdet(1:length(AllCellNames)));
tempMat(length(AllCellNames)+gap+1:end,1:length(AllCellNames)) = ConnMatrixPre(matIndex_commdet(length(AllCellNames)+1:end),matIndex_commdet(1:length(AllCellNames)));
tempMat(length(AllCellNames)+gap+1:end,1:length(AllCellNames)) = ConnMatrixPre(matIndex_commdet(length(AllCellNames)+1:end),matIndex_commdet(1:length(AllCellNames)));

%tempMat(1:length(AllCellNames),length(AllCellNames)+gap+1:end) = ConnMatrixPre(matIndex_commdet(1:length(AllCellNames)),matIndex_commdet(length(AllCellNames)+1:end));
%tempMat(length(AllCellNames)+gap+1:end,length(AllCellNames)+gap+1:end) = ConnMatrixPre(matIndex_commdet(length(AllCellNames)+1:end),matIndex_commdet(length(AllCellNames)+1:end));

figure;
%subplot(2,2,1)
%cspy(ConnMatrixPre(matIndex_commdet,matIndex_commdet),'Colormap',colorcet('R3'),'Levels',255,'MarkerSize',12);
tempMatThreshold = tempMat;
tempMatThreshold(tempMatThreshold>=5) =5;
cspy(tempMatThreshold,'Colormap',colorcet('R3','N',5),'Levels',5,'MarkerSize',8);
c = colorbar;
 c.Ticks = [1:5];
 c.TickLabels = [1:5];
hold on;

for i = 1:length(blocks)
    line([0,size(connMat_syn,2)]+0.5,[blocks(i),blocks(i)]+0.5,'color','k');
    line([blocks(i),blocks(i)]+0.5,[0,size(connMat_syn,2)]+0.5,'color','k');
end

line([0,size(connMat_syn,2)],[size(connMat_syn,2),size(connMat_syn,2)],'color','k');
line([size(connMat_syn,2),size(connMat_syn,2)],[0,size(connMat_syn,2)],'color','k');
line([0,size(connMat_syn,2)],[size(connMat_syn,2),size(connMat_syn,2)]+gap,'color','k');

set(gca,'XTick',[],'YTick',[]);
box on;

locateRS = find(ismember(AllCellNames(idx),RSall));
text(locateRS,repmat(-20,1,numel(locateRS)), sprintf('\\downarrow'), 'HorizontalAlignment','center', 'FontWeight','bold','Color','k');

