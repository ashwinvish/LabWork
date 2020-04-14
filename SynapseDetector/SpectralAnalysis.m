% spectral analysis
clear;
addpath(genpath('/Users/ashwin/Documents/'));
startup;colorSchemes

load ConnMatrixPre.mat
load ConnMatSummed.mat
load AllCells.mat
load cellIDs.mat
LoadDataFrame

confirmedALX = [76181 76201 76187 76184 76192 76197];
confirmedDBX = [76199 76200 76182 76183 76185 76186 76189 76191 76188 ];
confirmedBARHL = [76198 76190 76193 76194 76195 76196 ];
confirmedIntegrators = [confirmedALX,confirmedDBX,confirmedBARHL];
IO = [77483,77484,77485,77486,77487,77488,77489,77490,77491,77492,77493,77494,...
    77495,77496,77497,76353,77499,77500,77501,77502,77503,77504,77505,77506,77507,...
    77508,77509,77510,77511,77512,77513,77514,77515,77516,77517,77518,77519,77520,...
    77521,77522,77523,77524,77525,77526,77527,77528,77529,77530,77531,77532,77533,...
    77534,77535,77536,77537,77538,77539,78420,76350,76271,76448,76487,76355,76351,...
    76273,76272,76488,76498,76261,76266,78395,78401];

IO = IO(ismember(IO,AllCells));
%%
fileName = 'ranking_original_noMirror.csv';
fid = fopen(fileName);
% original pos,
% cellID,fanOut,fanIn,toIntegrator,fromInteg,toMotor,class,centrality
SpectralWithRS.complete = textscan(fid,'%d %d %d %d %d %d %d %s %f','Delimiter',',');
fclose(fid);

% [a,b] = sort(SpectralWithRS.complete{3}.*SpectralWithRS.complete{4},'descend'); % degree
% 
%  [m,n,v] = find(ConnMatrixPre);                    % m -postSynapse,row, n - column, v - number of synapses,value
%   G = digraph(n',m',v',size(AllCells,1)); % digraph(source, target, weight, numNodes)
%   G.Nodes.Names = num2str(AllCells);      % nodeNames
%   C = centrality(G,'betweenness');
%  [~,Csort] = sort(C,'descend');
%  clear G

%% plot without mirroring


ranksToConsider = 1:300; % eigenvector centrality
%ranksToConsider_degree = b(1:300); % degree
%[~,ranksToConsider] = ismember(AllCells(Csort(1:300)),SpectralWithRS.complete{2}); % betweenness

cellID  = SpectralWithRS.complete{2}(ranksToConsider);

findr456I = strcmp(SpectralWithRS.complete{8}(ranksToConsider),'Integrator_r456I');
findr456M = strcmp(SpectralWithRS.complete{8}(ranksToConsider),'Integrator_r456M');
findr456MI = strcmp(SpectralWithRS.complete{8}(ranksToConsider),'Integrator_r456MI');
findr56vm = strcmp(SpectralWithRS.complete{8}(ranksToConsider),'Integrator_r56M');
findr78ipsi = strcmp(SpectralWithRS.complete{8}(ranksToConsider),'Integrator_r78ipis');
findr78contra = strcmp(SpectralWithRS.complete{8}(ranksToConsider),'Integrator_r78contra');
findDO = strcmp(SpectralWithRS.complete{8}(ranksToConsider),'Vestibular_DO');
findRS = strcmp(SpectralWithRS.complete{8}(ranksToConsider),'RS');
findRS(find(cellID == 77099)) = 1; % MiT
findRS(find(cellID == 76202)) = 1; %M

contras = [81362,81599,81777,79400,79372,81337,78841,81671,81434,78983,78665];
%findContra = strcmp(SpectralWithRS.complete{8}(1:ranksToConsider),'Contra');
findContra = ismember(cellID,contras);
orphans = [77442,77120,78087,79394,79304,80211,81614,79242,77452,81608,78583,81605];
%findOrphans = strcmp(SpectralWithRS.complete{8}(1:ranksToConsider),'orphan');
findOrphans = ismember(cellID,orphans);
%findipsi_contra = strcmp(SpectralWithRS.complete{8}(1:ranksToConsider),'ipsi/');
ipsi_contra = [79510,76210,76543,76468,77927,76551,76614,79170,79584,77101,81431,79582,79362,79355,...
    78104,79214,79367,79356,79526,77232,79341,79378];
findipsi_contra = ismember(cellID,ipsi_contra);

 allIDed = logical(findr456I+findr456M+findr456MI+findr56vm+findr78ipsi+findr78contra+findDO+findRS+findContra+findOrphans+findipsi_contra);
 allUnIDed = logical(1-allIDed);
 cellID_cleaned = [cellID(findDO);cellID(findr456I);cellID(findr456M);cellID(findr456MI);cellID(findr56vm);cellID(findr78ipsi);cellID(findr78contra);cellID(findipsi_contra);cellID(findRS);cellID(allUnIDed)];

cellID_cleaned_noRS = [cellID(findDO);cellID(findr456I);cellID(findr456M);cellID(findr456MI);cellID(findr56vm);cellID(findr78ipsi);cellID(findr78contra);cellID(findipsi_contra);cellID(allUnIDed)];


% load ABDPutativeSaccadic.mat
% load ABDiPutativeSaccadic.mat
% 
% % IBNall = [82714,82710,79083,82264,78550,79084,82712,82709,82713,83184,77125,77942,77231,77247,77153,78685,77137,83183];
% DOs = [76631,76700,76656,76655,76660,76679,77393,77395,77375,...
%      77815,77803,77807,80591,80579,77255,77697,77364,80223,81126,80760,81117,...
%      80672,80622,80572,80569,81122,78426,81328,81575,81582,81870,82161,82165,81553,82169];
% % MOs = [76664 76784 76783 77394 77326 77772 77766 77756 77753 77767 77775 ...
% %      77780 80883 80247 80669 80698 80762 80748 80696 80761 80749 81070 81044 ...
% %      80991 80967 81008 80883 81045 80986 81023 78647 78619 80247 81141];
%   
% ABDr_CellIDs = [82140,82145,82143,77648,82146,77302,77301,82194,77705,77710,77305,77672,77300,77709,77661];
% ABDc_CellIDs = [82213,77654,77646,82212,77295,81172,77657,82197,77682,77154, 77652,77658,77628,77292,77688,82195,77296];
% ABDIr_CellIDs = [78553,77150,77886,78547,77631,77158,77665,78552,77668,77618,77634,78556];
% ABDIc_CellIDs = [78574,79051,77148,77625,77144,77692,77641,79066,77640,77643];
% 
%  
% MatOrder = [DOs';...
%     cellID(findr456I);...
%     cellID(findr456M);...
%     cellID(findr456MI);...
%     cellID(findr56vm);...
%     cellID(findr78ipsi);...
%     cellID(findr78contra);...
%     cellID(allUnIDed);...
%     ABDr_CellIDs';...
%     ABDc_CellIDs';...
%     ABDIr_CellIDs';...
%     ABDIc_CellIDs'];
% 
% [~,MatIndex]= ismember(MatOrder,AllCells);   
% r78Int = cellID(findr78contra);
% allIntegrators = [cellID(findr456I);cellID(findr456M);cellID(findr456MI);cellID(findr56vm);cellID(findr78ipsi);cellID(findr78contra)];
% allIntegrators_OSI = OSI(allIntegrators,df);
% allMonoIntegrators = [allIntegrators(allIntegrators_OSI == 1);allIntegrators(allIntegrators_OSI == -1)];
% 
% PutInts = cellID(allUnIDed);
% 
% findEnds = [find(MatOrder == ABDiPutativeSaccadic.cellIDs(end));...
%             find(MatOrder == MOs(end));...
%             find(MatOrder == r78Int(end));...
%             find(MatOrder == PutInts(end))];
%                    
% 
% subplot(2,2,3)
% cspy(ConnMatrixPre(MatIndex,MatIndex),'Colormap',colorcet('R3','N',15),'Levels',7,'MarkerSize',7);
% lineEdges = [zeros(length(findEnds),1),repmat(length(MatOrder),length(findEnds),1),findEnds,findEnds];
% for i = 1:size(lineEdges,2)
% line([lineEdges(i,1),lineEdges(i,2)],[lineEdges(i,3),lineEdges(i,4)],'color','k');
% line([lineEdges(i,3),lineEdges(i,4)],[lineEdges(i,1),lineEdges(i,2)],'color','k');
% end
% box on;
% title(sprintf('%s',fileName));
% 

%%

index = find(ismember(AllCells, cellID_cleaned ));
%connMat = ConnMatSummed(index,index);
connMat_syn = ConnMatrixPre(index,index);
connMat_sum = ConnMatSummed(index,index);
[m,n,v] = find(connMat_syn);                    % m -postSynapse,row, n - column, v - number of synapses,value
AllCellNames = AllCells(index);
% 
G = digraph(n',m',v',size(AllCellNames,1)); % digraph(source, target, weight, numNodes)
G.Nodes.Names = num2str(AllCellNames);      % nodeNames
nn = numnodes(G);
[s,t] = findedge(G);
A = sparse(s,t,G.Edges.Weight,nn,nn);


% gamma = 0.1:0.05:2;
% meanMod = [];
% meanClust = [];
% Modu = [];
% meanRandMod = [];
     jj = 1
 for kk = 1:500
%      randMat = randomize_graph_partial_und(A,randi([0,1],size(A)),1000); %
%      for jj = 1:length(gamma)
        Q0 = -1; Q1 = 0;                            % initialize modularity values
        i =1 ;
        while abs(Q1-Q0)>1e-5                        % while modularity increases
            Q0 = Q1;                                  % perform community detection
            [M, Q1] = community_louvain(A,0.38,1:size(connMat_syn,1));
            %[M2,Q2] = community_louvain(randMat,gamma(jj),1:size(randMat,1));
            i = i+1;
        end
%         [kk,jj]
%         
%                 P = participation_coef(connMat_syn,M);
%                   Z =  module_degree_zscore(connMat_syn,M,3);
%         %
%         %         %         symmWtMat = connMat_syn +connMat_syn';
%         %         %         m = sum(sum(symmWtMat));
%         %         %         MOD = 0;
%         %         %         COMu = unique(M);
%         %         %         for j=1:length(COMu)
%         %         %             Cj = find(M==COMu(j));
%         %         %             Ec = sum(sum(symmWtMat(Cj,Cj)));
%         %         %             Et = sum(sum(symmWtMat(Cj,:)));
%         %         %             MOD = MOD + Ec/m-(Et/m)^2;
%         %         %         end
%         % %         %Modu(kk,jj) = MOD;
         meanMod(kk,jj) = Q1;
         cluster(:,kk) = M;
         meanClust(kk,jj) = length(unique(M));
        % meanRandMod(kk,jj) = Q2;
%         % %
     end
 %end
%end
% % % % 
% subplot(4,4,1)
% errorbar(0.1:0.1:2,mean(meanMod),std(Modu),'-');
% hold on;
% scatter(0.1:0.1:2,mean(meanMod),5*mean(meanClust));
% xlabel('gamma');
% ylabel('Modularity');
% axis square
% box off;
% 
% subplot(4,4,2)
% scatter(mean(meanMod),mean(Modu))
% xlabel('modified modularity');
% ylabel('classic modularity');
% axis square
% box off;



%[a,b] = cluster_jl(connMat_syn,1,0);

%         symmWtMat = (connMat_syn + connMat_syn.')/2;
%         m = sum(sum(symmWtMat));
%         MOD = 0;
%         COMu = unique(a.COM{1});
%         newMat = zeros(COMu);
%         for j=1:length(COMu)
%             Cj = find(a.COM{1}==COMu(j));
%             Ec = sum(sum(symmWtMat(Cj,Cj)));
%             Et = sum(sum(symmWtMat(Cj,:)));
%             newMat(j,j) = sum(sum(symmWtMat(Cj,Cj)));
%             
%             for k = 1:length(COMu)
%                 Ck = find(a.COM{1}==COMu(k));
%                 newMat(j,k) = sum(sum(symmWtMat(Cj,Ck)));
%             end
%             MOD = MOD + Ec/m-(Et/m)^2;
%         end


% [X,Y,INDSORT] = grid_communities(M); % call function
% imagesc(connMat_syn(INDSORT,INDSORT));           % plot ordered adjacency matrix
%  hold on;                                 % hold on to overlay community visualization
%  plot(X,Y,'r','linewidth',2);             % plot community boundaries
%     
% association matrix
D = agreement(cluster); 
% random association matrix;
Tr = randi([min(unique(cluster)),max(unique(cluster))],[273,kk]);
Dr = agreement(Tr);
% thresholded agreement matrix;
T = consensus_und(D,max(Dr(:)),100);


%% figure 1 2 blocks only gamma  = 0.3
 [Vsort,idx] = sort(M);
 
% 
% block1.cellIDs = AllCellNames(idx(Vsort==1));
% block2.cellIDs = AllCellNames(idx(Vsort==2));
% 

%[Vsort,idx] = sort(a.COM{1});
% 
%  block1.cellIDs = AllCellNames(idx(Vsort==2));
%  block2.cellIDs = AllCellNames(idx(Vsort==1));
% 



%subplot(2,2,1)
tempMat = connMat_syn(idx,idx);
tempMapThreshold = tempMat;
%figure;
% subplot(4,4,1)
% histogram(tempMapThreshold,1:22,'EdgeColor','k','LineWidth',2,'FaceColor','w');
% box off;
% axis square;
% set(gca,'YScale','log');
% xlabel('Synapse count');
% ylabel('Number of connections');


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

locateRS = find(ismember(AllCellNames(idx),cellID(findRS)));
text(locateRS,repmat(-20,1,numel(locateRS)), sprintf('\\downarrow'), 'HorizontalAlignment','center', 'FontWeight','bold','Color','k');

% locateInts = find(ismember(AllCellNames(idx),confirmedIntegrators));
% text(locateInts,repmat(-10,1,numel(locateInts)), sprintf('\\downarrow'), 'HorizontalAlignment','center', 'FontWeight','bold','Color','k');
%% Add motor neurons
ABDr_CellIDs = [82140,82145,82143,77648,82146,77302,77301,82194,77705,77710,77305,77672,77300,77709,77661];
ABDc_CellIDs = [82213,77654,77646,82212,77295,81172,77657,82197,77682,77154, 77652,77658,77628,77292,77688,82195,77296];
ABDIr_CellIDs = [78553,77150,77886,78547,77631,77158,77665,78552,77668,77618,77634,78556];
ABDIc_CellIDs = [78574,79051,77148,77625,77144,77692,77641,79066,77640,77643];

RScells = cellID(findRS);
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


% for i = 1:length(blocks)
%     line([0,size(tempMat,2)]+0.5,[blocks(i),blocks(i)]+0.5,'color','k');
%     line([blocks(i),blocks(i)]+0.5,[0,size(tempMat,2)]+0.5,'color','k');
% end

%  line([0,size(tempMat,2)]+0.5,[size(Vsort,1),size(Vsort,1)]+0.5,'color','k');
%  line([size(Vsort,1),size(Vsort,1)]+0.5,[0,size(tempMat,2)]+0.5,'color','k');
%  
%  line([0,size(tempMat,2)]+0.5,[size(Vsort,1)+gap,size(Vsort,1)+gap]+0.5,'color','k');
%  line([size(Vsort,1)+gap,size(Vsort,1)+gap]+0.5,[0,size(tempMat,2)]+0.5,'color','k');
% % 
%  line([0,size(tempMat,2)+0.5],[size(Vsort,1)+32+gap,size(Vsort,1)+32+gap]+0.5,'color','k');
%  line([size(Vsort,1)+32+gap,size(Vsort,1)+32+gap]+0.5,[0,size(tempMat,2)+0.5],'color','k');


set(gca,'XTick',[],'YTick',[]);
box on;

% locateInts = find(ismember(AllCellNames(idx),confirmedIntegrators));
% text(locateInts,repmat(-10,1,numel(locateInts)), sprintf('\\downarrow'), 'HorizontalAlignment','center', 'FontWeight','bold','Color','k');

locateRS = find(ismember(AllCellNames(idx),cellID(findRS)));
text(locateRS,repmat(-20,1,numel(locateRS)), sprintf('\\downarrow'), 'HorizontalAlignment','center', 'FontWeight','bold','Color','k');

%% sub cluster only oculomotor block;

oculoMotorblock = AllCellNames(idx(Vsort==2));

index = find(ismember(AllCells, oculoMotorblock));
%connMat = ConnMatSummed(index,index);
oculoMotor_syn = ConnMatrixPre(index,index);
clear m;
clear n;
clear v;

[m,n,v] = find(oculoMotor_syn);                    % m -postSynapse,row, n - column, v - number of synapses,value
oculoMotorCellIDs = AllCells(index);
% 

G = digraph(n',m',v',size(oculoMotorCellIDs,1)); % digraph(source, target, weight, numNodes)
G.Nodes.Names = num2str(oculoMotorCellIDs);      % nodeNames
nn = numnodes(G);
[s,t] = findedge(G);
A = sparse(s,t,G.Edges.Weight,nn,nn);

% choose the correct gamma

%  gamma = 0.1:0.05:2;
%  meanMod = [];
%  meanRandMod = [];
% 
% for kk = 1:20
%     randMat = randomize_graph_partial_und(A,randi([0,1],size(A)),1000); %
%     parfor jj = 1:length(gamma)
%         Q0 = -1; Q1 = 0;                            % initialize modularity values
%         i =1 ;
%         while abs(Q1-Q0)>1e-5                        % while modularity increases
%             Q0 = Q1;                                  % perform community detection
%             [M, Q1] = community_louvain(A,gamma(jj),1:size(oculoMotor_syn,1));
%             [M2,Q2] = community_louvain(randMat,gamma(jj),1:size(randMat,1));
%             i = i+1;
%         end
%         meanMod(kk,jj) = Q1;
%         meanRandMod(kk,jj) = Q2;
%     end
% end
%  
% figure; 
% shadedErrorBar(gamma,mean(meanMod-meanRandMod),std(meanMod-meanRandMod)); axis square;

clear cluster;
jj = 1
 for kk = 1:500

        Q0 = -1; Q1 = 0;                            % initialize modularity values
        i =1 ;
        while abs(Q1-Q0)>1e-5                        % while modularity increases
            Q0 = Q1;                                  % perform community detection
            [M, Q1] = community_louvain(A,0.75,1:size(oculoMotor_syn,1));
            %[M2,Q2] = community_louvain(randMat,gamma(jj),1:size(randMat,1));
            %                      plot(i,Q1,'o');
            %                      hold on;
            %                      drawnow;
            i = i+1;
        end
         meanMod(kk,jj) = Q1;
         cluster(:,kk) = M;
         meanClust(kk,jj) = length(unique(M));
 end
    
% association matrix
oculomotor_D = agreement(cluster); 
% random association matrix;
oculomotor_Tr = randi([min(unique(cluster)),max(unique(cluster))],[size(oculoMotor_syn,1),kk]);
oculomotor_Dr = agreement(oculomotor_Tr);
% thresholded agreement matrix;
oculomotor_T = consensus_und(oculomotor_D,max(oculomotor_Dr(:)),100);

%% plot submodules

ABDr_CellIDs = [82140,82145,82143,77648,82146,77302,77301,82194,77705,77710,77305,77672,77300,77709,77661];
ABDc_CellIDs = [82213,77654,77646,82212,77295,81172,77657,82197,77682,77154, 77652,77658,77628,77292,77688,82195,77296];
ABDIr_CellIDs = [78553,77150,77886,78547,77631,77158,77665,78552,77668,77618,77634,78556];
ABDIc_CellIDs = [78574,79051,77148,77625,77144,77692,77641,79066,77640,77643];


[Vsort,idx] = sort(oculomotor_T);

matOrder_commdet_subMod = [oculoMotorCellIDs(idx)',...
    ABDr_CellIDs,ABDc_CellIDs,ABDIr_CellIDs,ABDIc_CellIDs];
   
[~,matIndex_commdet_subMod] = ismember(matOrder_commdet_subMod,AllCells);
gap = 20;
tempMat = zeros(length(matOrder_commdet_subMod)+gap);

tempMat(1:length(oculoMotorCellIDs),1:length(oculoMotorCellIDs)) = ConnMatrixPre(matIndex_commdet_subMod(1:length(oculoMotorCellIDs)),matIndex_commdet_subMod(1:length(oculoMotorCellIDs)));
tempMat(length(oculoMotorCellIDs)+gap+1:end,1:length(oculoMotorCellIDs)) = ConnMatrixPre(matIndex_commdet_subMod(length(oculoMotorCellIDs)+1:end),matIndex_commdet_subMod(1:length(oculoMotorCellIDs)));
%tempMat(1:length(AllCellNames),length(AllCellNames)+gap+1:end) = ConnMatrixPre(matIndex_commdet(1:length(AllCellNames)),matIndex_commdet(length(AllCellNames)+1:end));
%tempMat(length(AllCellNames)+gap+1:end,length(AllCellNames)+gap+1:end) = ConnMatrixPre(matIndex_commdet(length(AllCellNames)+1:end),matIndex_commdet(length(AllCellNames)+1:end));
clear tempMatThreshold;
tempMatThreshold = tempMat;
tempMatThreshold(tempMatThreshold>=5) =5;

figure;
cspy(tempMatThreshold,'Colormap',colorcet('R3','N',5),'Levels',5,'MarkerSize',8);
c = colorbar;
 c.Ticks = [1:5];
 c.TickLabels = [1:5];
hold on;

line([0,size(oculoMotor_syn,2)],[size(oculoMotor_syn,2),size(oculoMotor_syn,2)],'color','k');
line([size(oculoMotor_syn,2),size(oculoMotor_syn,2)],[0,size(oculoMotor_syn,2)],'color','k');
line([0,size(oculoMotor_syn,2)],[size(oculoMotor_syn,2),size(oculoMotor_syn,2)]+gap,'color','k');

line([57,57],[0,size(oculoMotor_syn,2)],'color','k');
line([0,size(oculoMotor_syn,2)],[57,57],'color','k');

line([57,57],[gap+size(oculoMotor_syn,2),size(tempMat,2)],'color','k');
line([0,size(oculoMotor_syn,2)],[size(oculoMotor_syn,2),size(oculoMotor_syn,2)]+gap+32,'color','k');


set(gca,'XTick',[],'YTick',[]);
box on;


%% Organize by centrality
% block1.centrality = SpectralWithRS.complete{9}(ismember(SpectralWithRS.complete{2},block1.cellIDs));
% block2.centrality = SpectralWithRS.complete{9}(ismember(SpectralWithRS.complete{2},block2.cellIDs));
% block3.centrality = SpectralWithRS.complete{9}(ismember(SpectralWithRS.complete{2},block3.cellIDs));


%% Submodules Axial blocks

AxialBlock = AllCellNames(idx(Vsort==1));

index = find(ismember(AllCells, AxialBlock));
%connMat = ConnMatSummed(index,index);
Axial_syn = ConnMatrixPre(index,index);
clear m;
clear n;
clear v;

[m,n,v] = find(Axial_syn);                    % m -postSynapse,row, n - column, v - number of synapses,value
AxialCellIDs = AllCells(index);
% 

G = digraph(n',m',v',size(AxialCellIDs,1)); % digraph(source, target, weight, numNodes)
G.Nodes.Names = num2str(AxialCellIDs);      % nodeNames
nn = numnodes(G);
[s,t] = findedge(G);
A = sparse(s,t,G.Edges.Weight,nn,nn);

% choose the correct gamma

 gamma = 0.1:0.05:2;
 meanMod = [];
 meanRandMod = [];

for kk = 1:20
    randMat = randomize_graph_partial_und(A,randi([0,1],size(A)),1000); %
    parfor jj = 1:length(gamma)
        Q0 = -1; Q1 = 0;                            % initialize modularity values
        i =1 ;
        while abs(Q1-Q0)>1e-5                        % while modularity increases
            Q0 = Q1;                                  % perform community detection
            [M, Q1] = community_louvain(A,gamma(jj),1:size(Axial_syn,1));
            [M2,Q2] = community_louvain(randMat,gamma(jj),1:size(randMat,1));
            i = i+1;
        end
        meanMod(kk,jj) = Q1;
        meanRandMod(kk,jj) = Q2;
    end
end
 
figure; 
shadedErrorBar(gamma,mean(meanMod-meanRandMod),std(meanMod-meanRandMod)); axis square;
%%
clear cluster;
jj = 1
 for kk = 1:500
        Q0 = -1; Q1 = 0;                            % initialize modularity values
        i =1 ;
        while abs(Q1-Q0)>1e-5                        % while modularity increases
            Q0 = Q1;                                  % perform community detection
            [M, Q1] = community_louvain(A,0.7,1:size(Axial_syn,1)); % 0.7 optimal
            %[M2,Q2] = community_louvain(randMat,gamma(jj),1:size(randMat,1));
            %                      plot(i,Q1,'o');
            %                      hold on;
            %                      drawnow;
            i = i+1;
        end
         meanMod(kk,jj) = Q1;
         cluster(:,kk) = M;
         meanClust(kk,jj) = length(unique(M));
 end
    
% association matrix
Axial_D = agreement(cluster); 
% random association matrix;
Axial_Tr = randi([min(unique(cluster)),max(unique(cluster))],[size(Axial_syn,1),kk]);
Axila_Dr = agreement(Axial_Tr);
% thresholded agreement matrix;
Axial_T = consensus_und(Axial_D,max(Axila_Dr(:)),100);


%%

[Vsort,idx] = sort(Axial_T);
matOrder_commdet_AxialsubMod = [AxialCellIDs(idx)'];
[~,matIndex_commdet_AxialsubMod] = ismember(matOrder_commdet_AxialsubMod,AllCells);

clear tempMat;


tempMat(1:length(AxialCellIDs),1:length(AxialCellIDs)) = ConnMatrixPre(matIndex_commdet_AxialsubMod(1:length(AxialCellIDs)),matIndex_commdet_AxialsubMod(1:length(AxialCellIDs)));
clear tempMatThreshold;
tempMatThreshold = tempMat;
tempMatThreshold(tempMatThreshold>=5) =5;

RScells = cellID(findRS);

figure;
cspy(tempMatThreshold,'Colormap',colorcet('R3','N',5),'Levels',5,'MarkerSize',12);
c = colorbar;
 c.Ticks = [1:5];
 c.TickLabels = [1:5];
hold on;

locs =find(diff(Vsort));

for i = 1: length(locs)
line([0,size(Axial_syn,2)],[locs(i),locs(i)],'color','k');
line([locs(i),locs(i)],[0,size(Axial_syn,2)],'color','k');
end
AxialCellIDs_sorted = AxialCellIDs(idx);
locateRS = find(ismember(AxialCellIDs_sorted,cellID(findRS)));

text(locateRS,repmat(-10,1,numel(locateRS)), sprintf('\\downarrow'), 'HorizontalAlignment','center', 'FontWeight','bold','Color','k');
box on;
%% blue, black, red
temp = distinguishable_colors(10,'w');

transform_swc_AV([77099,77265,77441,77449,77456,79961,80327,81410,82218],temp(5,:),[],true,false);
hold on;
transform_swc_AV([76202,79395,82217],temp(6,:),[],false,false);
transform_swc_AV([76562,77695,77931,78566,78577,79244],temp(7,:),[],false,false);



%%
temp = distinguishable_colors(10,'w');
col2 = [0,0.5,1]; % purples
col1 = [1,0.5,0] ; % Mustard

% hex2rgb('#9ec8cf') % bluish green
% hex2rgb('#7c6ecc') % purplish blue


figure(2) % only cell  example cells 76626,80629
transform_swc_AV(block1.cellIDs,col1,[],true,false);
hold on;
transform_swc_AV(block2.cellIDs,col2,[],false,false);

figure(3) % 
transform_swc_AV(block1.cellIDs,col1,[],true,false);
hold on;
transform_swc_AV(block2.cellIDs,col2,[],false,false);

% remove 76202 from block1;
block1.cellIDs(block1.cellIDs ==76202) = [];
block1.Origins = getOrigin(block1.cellIDs);
block2.Origins = getOrigin(block2.cellIDs);

figure(4)
Xvar = [block1.Origins(:,1);block2.Origins(:,1)];
Yvar = [block1.Origins(:,2);block2.Origins(:,2)];

grpIDs = [ones(size(block1.Origins,1),1);2*ones(size(block2.Origins,1),1)];
h = scatterhist(Xvar,Yvar,'Group',grpIDs,'color',[col1;col2],'Kernel','on');
set(h(1),'YDir','reverse','YLim',[400,800],'XLim',[250,400]);
set(h(3),'XDir','reverse');
daspect([1,1,1]);


% %%  with mirroring
% 
% fileName = 'ranking_original.csv';
% fid = fopen(fileName);
% SpectralWithRS.withMirroring = textscan(fid,'%d %d %d %d %d %d %d %s %f','Delimiter',',');
% fclose(fid);
% 
% ranksToConsider = 300;
% cellID  = SpectralWithRS.withMirroring{2}(1:ranksToConsider);
% 
% findr456I = strcmp(SpectralWithRS.withMirroring{8}(1:ranksToConsider),'Integrator_r456I');
% findr456M = strcmp(SpectralWithRS.withMirroring{8}(1:ranksToConsider),'Integrator_r456M');
% findr456MI = strcmp(SpectralWithRS.withMirroring{8}(1:ranksToConsider),'Integrator_r456MI');
% findr56vm = strcmp(SpectralWithRS.withMirroring{8}(1:ranksToConsider),'Integrator_r56M');
% findr78ipsi = strcmp(SpectralWithRS.withMirroring{8}(1:ranksToConsider),'Integrator_r78ipis');
% findr78contra = strcmp(SpectralWithRS.withMirroring{8}(1:ranksToConsider),'Integrator_r78contra');
% findDO = strcmp(SpectralWithRS.withMirroring{8}(1:ranksToConsider),'Vestibular_DO');
% findRS = strcmp(SpectralWithRS.withMirroring{8}(1:ranksToConsider),'RS');
% contras = [81362,81599,81777,79400,79372,81337,78841,81671,81434];
% %findContra = strcmp(SpectralWithRS.complete{8}(1:ranksToConsider),'Contra');
% findContra = ismember(cellID,contras);
% orphans = [77442,77120,78087,79394,79304,80211,81614,79242,77452,81608,78583,81605];
% %findOrphans = strcmp(SpectralWithRS.complete{8}(1:ranksToConsider),'orphan');
% findOrphans = ismember(cellID,orphans);
% %findipsi_contra = strcmp(SpectralWithRS.complete{8}(1:ranksToConsider),'ipsi/');
% ipsi_contra = [79510,76210,76543,76468,77927,76551,76614,79170,79584,77101,81431,79582,79362,79355,...
%     78104,79214,79367,79356,79526,77232,79341];
% findipsi_contra = ismember(cellID,ipsi_contra);
% 
% allIDed = logical(findr456I+findr456M+findr456MI+findr56vm+findr78ipsi+findr78contra+findDO+findRS+findContra+findOrphans+findipsi_contra);
% allUnIDed = logical(1-allIDed);
% cellID_cleaned = [cellID(findr456I);cellID(findr456M);cellID(findr456MI);cellID(findr56vm);cellID(findr78ipsi);cellID(findr78contra);cellID(findipsi_contra);cellID(findRS);cellID(allUnIDed)];
% 
% index = find(ismember(AllCells,cellID_cleaned));
% connMat = ConnMatrixPre(index,index);
% [m,n,v] = find(connMat);                    % m -postSynapse,row, n - column, v - number of synapses,value
% AllCellNames = AllCells(index);
% 
% G = digraph(n',m',v',size(AllCellNames,1)); % digraph(source, target, weight, numNodes)
% G.Nodes.Names = num2str(AllCellNames);      % nodeNames
% nn = numnodes(G);
% [s,t] = findedge(G);
% A = sparse(s,t,G.Edges.Weight,nn,nn);
% 
% for jj = 1:1
%     Q0 = -1; Q1 = 0;                            % initialize modularity values
%     i =1 ;
%     while abs(Q1-Q0)>1e-4                         % while modularity increases
%         Q0 = Q1;                                  % perform community detection
%         [M(:,jj), Q1(jj)] = community_louvain(A,0.4,1:size(A,1));
% %         plot(i,Q1-Q0,'o');
% %         hold on;
% %         drawnow;
%         i = i+1;
%     end
% end
% 
% [Vsort,idx] = sort(M);
% 
% block1.cellIDs = AllCellNames(idx(Vsort==1));
% block2.cellIDs = AllCellNames(idx(Vsort==2));
% 
% figure(1); 
% subplot(2,2,1)
% cspy(A(idx,idx),'Colormap',colorcet('L3'),'Levels',255,'MarkerSize',15);
% blocks = find(diff(round(Vsort)));
% c = colorbar;
% % c.Ticks = [0,0.2,0.4,0.6,0.8,1];
% % c.TickLabels = [0:0.01:0.02,0.03,0.04,0.05];
% hold on;
% 
% for i = 1:length(blocks)
%     line([0,size(A,2)]+0.5,[blocks(i),blocks(i)]+0.5,'color','k');
%     line([blocks(i),blocks(i)]+0.5,[0,size(A,2)]+0.5,'color','k');
% end
% box on;
% set(gca,'Xtick',[],'Ytick',[]);
% 


%%

figure; 
cspy(connMatForPlot_Norm,'Colormap',colorcet('R3'),'Levels',255,'MarkerSize',15);
blocks = find(diff(round(Vsort)));
c = colorbar;
c.Ticks = [0,0.2,0.4,0.6,0.8,1];
c.TickLabels = [0:0.01:0.02,0.03,0.04,0.05];
hold on;

for i = 1:length(blocks)
    line([0,size(matIndex_commdet,2)]+0.5,[blocks(i),blocks(i)]+0.5,'color','k');
    line([blocks(i),blocks(i)]+0.5,[0,size(matIndex_commdet,2)]+0.5,'color','k');
end
line([0,size(matIndex_commdet,2)],102+[size(cellID_cleaned,1),size(cellID_cleaned,1)],'color','k');
line(102+[size(cellID_cleaned,1),size(cellID_cleaned,1)],[0,size(matIndex_commdet,2)],'color','k');

line([0,size(matIndex_commdet,2)],[13,13],'color','k');
line([13,13],[0,size(matIndex_commdet,2)],'color','k');

line([0,size(matIndex_commdet,2)],[23,23],'color','k');
line([23,23],[0,size(matIndex_commdet,2)],'color','k');

line([0,size(matIndex_commdet,2)],[33,33],'color','k');
line([33,33],[0,size(matIndex_commdet,2)],'color','k');

line([0,size(matIndex_commdet,2)],[68,68],'color','k');
line([68,68],[0,size(matIndex_commdet,2)],'color','k');

line([0,size(matIndex_commdet,2)],[102,102],'color','k');
line([102,102],[0,size(matIndex_commdet,2)],'color','k');

line([0,size(matIndex_commdet,2)],102+32+[size(cellID_cleaned,1),size(cellID_cleaned,1)],'color','k');
line(102+32+[size(cellID_cleaned,1),size(cellID_cleaned,1)],[0,size(matIndex_commdet,2)],'color','k');


%set(gca,'color',[0.8,0.8,0.8],'XTick',[],'YTick',[]);
box on;

locateRS = find(ismember(AllCellNames(idx),cellID(findRS)));
locateMonoIntegrators = find(ismember(AllCellNames(idx),allMonoIntegrators));
locateBinoIntegrators = find(ismember(AllCellNames(idx),setdiff(allIntegrators,allMonoIntegrators)));
text(102+locateRS,repmat(-10,1,numel(locateRS)), sprintf('\\downarrow'), 'HorizontalAlignment','center', 'FontWeight','bold','Color','k');
text(repmat(-10,1,numel(locateRS)),102+locateRS, sprintf('\\rightarrow'), 'HorizontalAlignment','center', 'FontWeight','bold','Color','k');

temp = distinguishable_colors(10,'w');

figure;
startup
subplot(4,4,1)
histogram(OSI(Block3',df),-1:0.1:1,'DisplayStyle','stairs','EdgeColor',temp(7,:),'LineWidth',2);
axis square;
ylabel('count');
xlabel('OSI');
offsetAxes(gca);
box off;


% text(102+locateMonoIntegrators,repmat(-10,1,numel(locateMonoIntegrators)), sprintf('\\downarrow'), 'HorizontalAlignment','center', 'FontWeight','bold','Color','r');
% text(102+locateBinoIntegrators,repmat(-10,1,numel(locateBinoIntegrators)), sprintf('\\downarrow'), 'HorizontalAlignment','center', 'FontWeight','bold','Color','g');

%% 
%blocindex = [find(Vsort ==2);find(Vsort ==3)];

for i = 1:length(idx)
    %blocColors(i,:) = temp(Vsort(blocindex(i)),:);
    blocColors(i,:) = temp(3+Vsort(i),:);
end

subplot(1,3,1)
%transform_swc_AV(AllCellNames(idx(blocindex)),blocColors,[],true,false);
transform_swc_AV(AllCellNames(idx(Vsort==1)),temp(5,:),[],true,false);
subplot(1,3,2)
transform_swc_AV(AllCellNames(idx(Vsort==2)),temp(6,:),[],true,false);
subplot(1,3,3)
transform_swc_AV(AllCellNames(idx(Vsort==3)),temp(8,:),[],true,false);

colormap(blocColors);
c = colorbar;
c.Ticks = [0.1,0.4,0.5,0.7,0.9];
c.TickLabels = [1,2,3,4,5];


RScellsbloc1 = find(ismember(AllCellNames(find(Vsort ==1)),cellID(findRS)));
RScellsbloc4 = find(ismember(AllCellNames(find(Vsort ==4)),cellID(findRS)));

subplot(1,3,2)
transform_swc_AV(AllCellNames(RScellsbloc1),temp(1,:),[],true,false); % block 1
transform_swc_AV(AllCellNames(RScellsbloc4),temp(4,:),[],false,false); % block 3


% subplot(1,3,3)
% transform_swc_AV(AllCellNames(idx(115:139)),temp(2,:),[],true,false); % block 3


%%

spectralColors = cbrewer('seq','Blues',length(cellID));
transform_swc_AV(cellID,spectralColors,[],true,false);
colormap(spectralColors);
 c = colorbar;
 c.Ticks = [0,1];
 c.TickLabels = {'High','low'};
 



%% plot locations


transform_swc_AV(cellID(findr456I),SaccABDicolor,[],true,false);
hold on;

transform_swc_AV(cellID(findr456M),SaccABDcolor,[],false,false);
transform_swc_AV(cellID(findr456MI),[0.5,0.5,0.5],[],false,false);
transform_swc_AV(cellID(findr56vm),[0,0,1],[],true,false);
transform_swc_AV(cellID(findr78ipsi),ALXcolor,[],true,false);
transform_swc_AV(cellID(findr78contra),DBXcolor,[],false,false);

% transform_swc_AV(cellID(findDO),Vestcolor,[],true,false);
% transform_swc_AV(cellID(findRS),[1,0,0],[],true,false);

%spectralColors = cmap('R3','N',length(cellID(findPutative)));
spectralColors = distinguishable_colors(length(cellID(findPutative)),'w');
transform_swc_AV(cellID(findPutative),spectralColors,[],true,false);

%%
subplot(4,4,1)
scatter(SpectralWithRS.complete{9}(1:300),SpectralWithRS.complete{3}(1:300),5,'.','MarkerEdgeColor','k');
set(gca,'XScale','log');
ylabel('fan out')
xlabel('Centrality')
hold on
yyaxis right
scatter(SpectralWithRS.complete{9}(1:300),SpectralWithRS.complete{4}(1:300),5,'.','MarkerEdgeColor','r');
ylabel('fan in')
set(gca,'YColor','r');
axis square;


subplot(4,4,2)
cols = cmap('I3','N',300);
scatter(SpectralWithRS.complete{3}(1:300),SpectralWithRS.complete{4}(1:300),10,cols,'filled','MarkerEdgeColor','none');
colormap(cols)
c = colorbar;
c.Ticks = [0,1];
c.TickLabels = {'Low','High'};
c.Label.String = 'Centrality';
axis square;
f = showfit(ezfit(single(SpectralWithRS.complete{3}(1:300)),single(SpectralWithRS.complete{4}(1:300)),'affine'),'fitcolor','k','dispeqboxmode','off');
text(300,600,sprintf('r = %0.2f',f.r));
xlabel('fan in');
ylabel('fan out');

%%

 
TOmirror = [76257,77743,80741,81569,80738];
TOmirror = TOmirror(logical(isExistReRoot(TOmirror)));

Vest.cellIDs = [DOs';MOs';TOmirror'];
EBN.cellIDs = [ABDPutativeSaccadic.cellIDs';ABDiPutativeSaccadic.cellIDs'];

ABD_M = vertcat(ABDr_CellIDs', ABDc_CellIDs');
ABD_I = vertcat(ABDIr_CellIDs', ABDIc_CellIDs');

 % IBN mirror

IBN.mirrorPop = [77303,79799,80194,80847,80701,77355,77353,80737,80816,80707,81169,77855,79799,80700,80868];
IBN.mirrorPop  = IBN.mirrorPop(logical(isExistReRoot(IBN.mirrorPop)));

 
% % V/S ratio
allIntegrators.cellIDs = [PutInts;r78Int];
allIntegrators.cellIDs = allIntegrators.cellIDs(isExistReRoot(allIntegrators.cellIDs));
allIntegrators.Orgins = getOrigin(allIntegrators.cellIDs);
% 
for i = 1:size(allIntegrators.cellIDs,1)
    a = SynapticPartners(allIntegrators.cellIDs(i),1,df);
    vestSum(i) = sum(ismember(a,Vest.cellIDs));
    saccSum(i) = sum(ismember(a,EBN.cellIDs));
    clear a;
end

[~,sortorder] = sort(allIntegrators.Orgins(:,2));
cols = colorcet('D1A','N', size(allIntegrators.cellIDs,1));
VSratio = (vestSum-saccSum)./(vestSum+saccSum);
% 
figure
subplot(4,4,1)
scatter(allIntegrators.Orgins(:,2),VSratio);
view([90,90]);
axis square;
xlabel('C<-->R');
ylabel('V-S/V+S');

subplot(4,4,2)
scatter(vestSum(sortorder),saccSum(sortorder),50,cols(sortorder,:),'filled','MarkerEdgeColor','k');
axis square;
colormap(cols)
c = colorbar;
c.Ticks = [0,1];
c.TickLabels = {'Rostral','Caudal'};
xlabel('Vest. synapses')
ylabel('EBN synapses');
offsetAxes(gca);

% 


%% Assemble connectome for modeling

load ConnMatrixPre.mat
load AllCells.mat
load ABDPutativeSaccadic.mat
load ABDiPutativeSaccadic.mat

IBNall = [82714,82710,79083,82264,78550,79084,82712,82709,82713,83184,77125,77942,77231,77247,77153,78685,77137,83183];
DOs = [76631,76700,76656,76655,76660,76679,77393,77395,77375,...
     77815,77803,77807,80591,80579,77255,77697,77364,80223,81126,80760,81117,...
     80672,80622,80572,80569,81122,78426,81328,81575,81582,81870,82161,82165,81553,82169];
MOs = [76664 76784 76783 77394 77326 77772 77766 77756 77753 77767 77775 ...
     77780 80883 80247 80669 80698 80762 80748 80696 80761 80749 81070 81044 ...
     80991 80967 81008 80883 81045 80986 81023 78647 78619 80247 81141];
RS = SpectralWithRS.RS;
% Integrators_ipsi = SpectralWithRS.allIntegrators;
% Integrators_contra = [];
% Integrators_putative = [];

[~,b] = sort(SpectralWithRS.IntegratorsOSI);
Integrators_ipsi = SpectralWithRS.Integrators(b);
Integrators_contra = SpectralWithRS.r78Contra;
Integrators_putative = [SpectralWithRS.notIntegrators;SpectralWithRS.ipsiContra];

ABDr_CellIDs = [82140,82145,82143,77648,82146,77302,77301,82194,77705,77710,77305,77672,77300,77709,77661];
ABDc_CellIDs = [82213,77654,77646,82212,77295,81172,77657,82197,77682,77154, 77652,77658,77628,77292,77688,82195,77296];
ABDIr_CellIDs = [78553,77150,77886,78547,77631,77158,77665,78552,77668,77618,77634,78556];
ABDIc_CellIDs = [78574,79051,77148,77625,77144,77692,77641,79066,77640,77643];

matOrder = [ABDPutativeSaccadic.cellIDs';...
    ABDiPutativeSaccadic.cellIDs';...
    IBNall';...
    DOs';...
    MOs';...
    RS;...
    Integrators_ipsi;...
   Integrators_contra;...
   Integrators_putative;...
    ABDr_CellIDs';...
    ABDc_CellIDs';...
    ABDIr_CellIDs';...
    ABDIc_CellIDs'];


matCellIDs_postSpectralAnalysis = [repmat("EBN_m",length(ABDPutativeSaccadic.cellIDs),1);...
    repmat("EBN_i",length(ABDiPutativeSaccadic.cellIDs),1);...
    repmat("IBN",length(IBNall),1);...
    repmat("DO",length(DOs),1);...
    repmat("MO",length(MOs),1);...
    repmat("RS",length(RS),1);...
    repmat("Integ_ipsi",length(Integrators_ipsi),1);...
    repmat("Integ_contra",length(Integrators_contra),1);...
   repmat("Integ_put",length(Integrators_putative),1);...
    repmat("ABDr",length(ABDr_CellIDs),1);...
    repmat("ABDc",length(ABDc_CellIDs),1);...
    repmat("ABDIr",length(ABDIr_CellIDs),1);...
    repmat("ABDIc",length(ABDIc_CellIDs),1)];

for i = 1:length(matOrder)
matTotalInputs_postSpectralAnalysis(i) = length(SynapticPartners(matOrder(i),1,df));
end

[~,matIndex]= ismember(matOrder,AllCells);   

connMat_forAlex_postSpectralAnalysis = ConnMatrixPre(matIndex,matIndex);
% set diagonal elements to 0
connMat_forAlex_postSpectralAnalysis(logical(eye(size(connMat_forAlex_postSpectralAnalysis,1)))) = 0; 

for i = 1:size(connMat_forAlex_postSpectralAnalysis,1)
connMat_forAlex_postSpectralAnalysis_Normalized(i,:) = connMat_forAlex_postSpectralAnalysis(i,:)./matTotalInputs_postSpectralAnalysis(i);
end
connMat_forAlex_postSpectralAnalysis_Normalized(connMat_forAlex_postSpectralAnalysis_Normalized>0.05) = 0.05;
%%
subplot(2,2,1)
cspy(connMat_forAlex_postSpectralAnalysis_Normalized,'Colormap',colorcet('R3','N',15),'Levels',15,'MarkerSize',7);

end_EBN_m = length(ABDPutativeSaccadic.cellIDs);
end_EBN_i = end_EBN_m + length(ABDiPutativeSaccadic.cellIDs);
end_IBN = length(IBNall)+ end_EBN_i;
end_DO = length(DOs)+ end_IBN;
end_MO = length(MOs)+ end_DO;
end_RS = length(RS) + end_MO;
end_Integ_ipsi  = length(Integrators_ipis) + end_RS;
end_Integ_contra = length(Integrators_contra) + end_Integ_ipsi;
end_Integ_putative = length(Integrators_putative) + end_Integ_contra;
end_ABDm  = length([[ABDr_CellIDs],[ABDc_CellIDs]]) + end_Integ_putative;

line([0,size(connMat_forAlex_postSpectralAnalysis,1)]+0.5,[end_EBN_m,end_EBN_m]+0.5); % EBN_m
line([0,size(connMat_forAlex_postSpectralAnalysis,1)]+0.5,[end_EBN_i,end_EBN_i]+0.5); % EBN_i
line([0,size(connMat_forAlex_postSpectralAnalysis,1)]+0.5,[end_IBN,end_IBN]+0.5); % IBN
line([0,size(connMat_forAlex_postSpectralAnalysis,1)]+0.5,[end_DO,end_DO]+0.5); % DOs
line([0,size(connMat_forAlex_postSpectralAnalysis,1)]+0.5,[end_MO,end_MO]+0.5); % MOs
line([0,size(connMat_forAlex_postSpectralAnalysis,1)]+0.5,[end_RS,end_RS]+0.5); % RS
line([0,size(connMat_forAlex_postSpectralAnalysis,1)]+0.5,[end_Integ_ipsi,end_Integ_ipsi]+0.5); % Integ_ipsi
line([0,size(connMat_forAlex_postSpectralAnalysis,1)]+0.5,[end_Integ_contra,end_Integ_contra]+0.5); % Integ_contrs
line([0,size(connMat_forAlex_postSpectralAnalysis,1)]+0.5,[end_Integ_putative,end_Integ_putative]+0.5); % Putative ipntegrators
line([0,size(connMat_forAlex_postSpectralAnalysis,1)]+0.5,[end_ABDm,end_ABDm]+0.5); % EBN_m


line([end_EBN_m,end_EBN_m]+0.5,[0,size(connMat_forAlex_postSpectralAnalysis,1)]+0.5); % EBN_m
line([end_EBN_i,end_EBN_i]+0.5,[0,size(connMat_forAlex_postSpectralAnalysis,1)]+0.5); % EBN_i
line([end_IBN,end_IBN]+0.5,[0,size(connMat_forAlex_postSpectralAnalysis,1)]+0.5); % IBN
line([end_DO,end_DO]+0.5,[0,size(connMat_forAlex_postSpectralAnalysis,1)]+0.5); % DOs
line([end_MO,end_MO]+0.5,[0,size(connMat_forAlex_postSpectralAnalysis,1)]+0.5); % MOs
line([end_RS,end_RS]+0.5,[0,size(connMat_forAlex_postSpectralAnalysis,1)]+0.5); % RS
line([end_Integ_ipsi,end_Integ_ipsi]+0.5,[0,size(connMat_forAlex_postSpectralAnalysis,1)]+0.5); % Integ_ipsi
line([end_Integ_contra,end_Integ_contra]+0.5,[0,size(connMat_forAlex_postSpectralAnalysis,1)]+0.5); % Integ_contrs
line([end_Integ_putative,end_Integ_putative]+0.5,[0,size(connMat_forAlex_postSpectralAnalysis,1)]+0.5); % Putative ipntegrators
line([end_ABDm,end_ABDm]+0.5,[0,size(connMat_forAlex_postSpectralAnalysis,1)]+0.5); % EBN_m

set(gca,'XTick',[],'YTick',[],'color',[0.9,0.9,0.9]);
box on;
colormap(colorcet('R3','N',15));
colorbar('Ticks',[0,0.5,1],'TickLabels', [0,max(connMat_forAlex_postSpectralAnalysis_Normalized(:))/2,max(connMat_forAlex_postSpectralAnalysis_Normalized(:))],...
    'TickDirection','out');



% 
% save('connMat_forAlex_postSpectralAnalysis.mat','connMat_forAlex_postSpectralAnalysis');
% save('matTotalInputs_postSpectralAnalysis.mat','matTotalInputs_postSpectralAnalysis');
% save('matCellIDs_postSpectralAnalysis.mat','matCellIDs_postSpectralAnalysis');
% 

