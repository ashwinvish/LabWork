%clear;
load('MelanieDBXCells.mat');
load('MelanieALXCells.mat');
load('MelanieDeconvCells.mat');
load('STA.mat');
load('STA_raw.mat');
%load('FiringRates_Complete.mat');

OriginalCellOrderDBX = [8,9,10,11,21,12,15,2,3];
DbxCells = [76182,76183,76185,76186,76188,76189,76191,76199,76200];

STA = 0 % 0 for SVD computed Firing rates; 1 for STA only.

if STA == 1
    t = [-2:0.05:7]; % time period in sec
    % Using the STAs as is
    
    FiringRates = [STAall(:,OriginalCellOrderDBX) , DBX_vglut_neg/100 , DBX_vglut/100, ALX_neg/100, ALX_pos/100]; % divide by 100 to convert back from %
    FiringDBX = [STAall(:,OriginalCellOrderDBX) , DBX_vglut_neg/100 , DBX_vglut/100,];
    
    DBXpop = [DBX_vglut_neg/100 , DBX_vglut/100]; % melanies data is reported as %
    ALXpop = [ALX_pos/100, ALX_neg/100]; % melanies data.
    CellOrder = [ones(1,size(DbxCells,2)), 2* ones(1,size(DBX_vglut_neg,2)), 3* ones(1,size(DBX_vglut,2)), ...
        4* ones(1,size(ALX_neg,2)), 5* ones(1,size(ALX_pos,2))];
    
    
    % consider only the first 5 components of the trace
    [A,B,C,D,E,F] = pca(FiringRates);
    sprintf('first %d components capture %f of the data',2, sum(E(1:5)))
    FiringRates = B(:,1:5)*A(:,1:5)'+ F;
    
    clear A;
    clear B;
    clear C;
    clear D;
    clear E;
    clear F;
    [A,B,C,D,E,F] = pca(FiringDBX);
    sprintf('first %d components capture %f of the data',5, sum(E(1:5)))
    FiringDBX = B(:,1:5)*A(:,1:5)'+ F;
    
else
    % Using the svd computed Firing Rates.
    
    t = 0:0.05:0.05*(180);
    t=t-2;
    
    %FiringSTA = [STA_raw(:,OriginalCellOrderDBX) , DBX_neg_dstaf/100 , DBX_pos_dstaf/100, ALX_neg_dstaf/100, ALX_pos_dstaf/100]; % divide by 100 to convert back from %
    FiringSTA = [STAall(:,OriginalCellOrderDBX) , DBX_vglut_neg/100 , DBX_vglut/100, ALX_neg/100, ALX_pos/100]; % divide by 100 to convert back from %

    normFiringSTA = FiringSTA ./ max(FiringSTA);
    FiringRates = svd_rates(normFiringSTA, t, 1.9,[3,3],-100);
    %FiringRates = svd_rates(FiringSTA, t, 1.9,[3,3],-100);

    
    CellOrder = [ones(1,size(DbxCells,2)), 2* ones(1,size(DBX_vglut_neg,2)), 3* ones(1,size(DBX_vglut,2)), ...
        4* ones(1,size(ALX_neg,2)), 5* ones(1,size(ALX_pos,2))];
    
    FiringDBX = [FiringRates(:,OriginalCellOrderDBX) , FiringRates(:,find(CellOrder==2)), FiringRates(:,find(CellOrder==3))];
    DBXpop = [ FiringRates(:,find(CellOrder==2)), FiringRates(:,find(CellOrder==3))]; % melanies data is reported as %
    ALXpop = [ FiringRates(:,find(CellOrder==4)), FiringRates(:,find(CellOrder==5))]; % melanies data.
    
    
    
    % consider only the first 5 components of the trace
    [A,B,C,D,E,F] = pca(FiringRates);
    sprintf('first %d components capture %f of the data',2, sum(E(1:2)))
    FiringRates = B(:,1:2)*A(:,1:2)'+ F;
    
    clear A;
    clear B;
    clear C;
    clear D;
    clear E;
    clear F;
    
    [A,B,C,D,E,F] = pca(FiringDBX);
    sprintf('first %d components capture %f of the data',2, sum(E(1:2)))
    FiringDBX = B(:,1:2)*A(:,1:2)'+ F;
    
    t = t(1:end-1);
    
end
%% plot raw data

figure(1)
subplot(4,4,1)
shadedErrorBar(t,mean(FiringRates(:,CellOrder==2),2),std(FiringRates(:,CellOrder==2),[],2),'lineProps',{'b'});
xlabel('Time (sec)');
ylabel('df/f');
axis square;

subplot(4,4,1)
shadedErrorBar(t,mean(FiringRates(:,CellOrder==3),2),std(FiringRates(:,CellOrder==3),[],2),'lineProps',{'r'});
title('blue:DBXneg, red: DBXpos');
xlabel('Time (sec)');
ylabel('df/f');
axis square;

subplot(4,4,3)
shadedErrorBar(t,mean(FiringRates(:,CellOrder==4),2),std(FiringRates(:,CellOrder==4),[],2),'lineProps',{'b'});
xlabel('Time (sec)');
ylabel('df/f');
axis square;

subplot(4,4,3)
shadedErrorBar(t,mean(FiringRates(:,CellOrder==5),2),std(FiringRates(:,CellOrder==5),[],2),'lineProps',{'r'});
title('blue:ALXneg, red: ALXpos');
xlabel('Time (sec)');
ylabel('df/f');
axis square;


%% clustering analysis

% only melanies data

CorrValsALL = corr(FiringRates);
CorrValsALL = CorrValsALL - eye(size(CorrValsALL));

eva = evalclusters(FiringRates','kmeans','silhouette','Distance','correlation','KList',[1:6]);

subplot(4,4,5)
plot(eva);
title('Number of optimal clusters');
axis square;
box off;


if size(find(eva.OptimalY ==1),1)< size(find(eva.OptimalY ==2))
    Index2 = find(eva.OptimalY ==1);
    Index1 = find(eva.OptimalY ==2);
    Index3 = find(eva.OptimalY ==3)
else
    Index2 = find(eva.OptimalY ==2);
    Index1 = find(eva.OptimalY ==1);
    Index3 = find(eva.OptimalY ==3)

end

subplot(4,4,6)
FiringOrdered = [FiringRates(:,Index2),FiringRates(:,Index1),FiringRates(:,Index3)];
imagesc(FiringOrdered');
colormap(colorcet('L16'));
colorbar;
set(gca,'XTick',[0,41,81,121,161],'XTickLabel',[-2,0,2,4,6]);
xlabel('Time(sec)');
box off;
ylabel('Neurons');
daspect([1,1,1]);

subplot(4,4,7)
shadedErrorBar(t,mean(FiringRates(:,Index2),2),std(FiringRates(:,Index2),[],2),'lineprops',{'r'});
hold on;
shadedErrorBar(t,mean(FiringRates(:,Index1),2),std(FiringRates(:,Index1),[],2),'lineprops',{'b'});
shadedErrorBar(t,mean(FiringRates(:,Index3),2),std(FiringRates(:,Index3),[],2),'lineprops',{'g'});

axis square;
xlabel('Time(sec)');
ylabel('Df/f');
box off;

subplot(4,4,8)
histogram(CellOrder(Index2),'FaceColor','r');
hold on;
histogram(CellOrder(Index1),'FaceColor','b');
histogram(CellOrder(Index3),'FaceColor','g');

xticks([1,2,3,4,5]);
xticklabels({'Int','DBX-','DBX+','ALX-','ALX+'});
xtickangle(45);
axis square;
box off;

clear leadIndex;
clear lagIndex;

%% Make plot with DBX only;

CorrValsDBX = corr(FiringDBX);
CorrValsDBX = CorrValsDBX - eye(size(CorrValsDBX));

evaDBX = evalclusters(FiringDBX','kmeans','silhouette','Distance','correlation','KList',[1:6]);

subplot(4,4,9)
plot(evaDBX);
title('Number of optimal clusters');
axis square;
box off;


if size(find(evaDBX.OptimalY ==1),1)> size(find(evaDBX.OptimalY ==2))
    Index1 = find(evaDBX.OptimalY ==1);
    Index2 = find(evaDBX.OptimalY ==2);
    Index3 = find(evaDBX.OptimalY ==3);
else
    Index1 = find(evaDBX.OptimalY ==2);
    Index2 = find(evaDBX.OptimalY ==1);
    Index3 = find(evaDBX.OptimalY ==3);
end

subplot(4,4,10)
FiringOrderedDBX = [FiringDBX(:,Index2),FiringDBX(:,Index1),FiringDBX(:,Index3)];
imagesc(FiringOrderedDBX');
colormap(colorcet('L16'));
colorbar;
set(gca,'XTick',[0,41,81,121,161],'XTickLabel',[-2,0,2,4,6]);
xlabel('Time(sec)');
box off;
ylabel('Neurons');
daspect([1,1,1]);


subplot(4,4,11)
shadedErrorBar(t,mean(FiringDBX(:,Index2),2),std(FiringDBX(:,Index2),[],2),'lineprops',{'r'});
hold on;
shadedErrorBar(t,mean(FiringDBX(:,Index1),2),std(FiringDBX(:,Index1),[],2),'lineprops',{'b'});
shadedErrorBar(t,mean(FiringDBX(:,Index3),2),std(FiringDBX(:,Index3),[],2),'lineprops',{'g'});

axis square;
xlabel('Time(sec)');
ylabel('Df/f');
box off;

subplot(4,4,12)
histogram(CellOrder(Index2),'FaceColor','r');
hold on;
histogram(CellOrder(Index1),'FaceColor','b');
histogram(CellOrder(Index3),'FaceColor','g');

xticks([1,2,3]);
xticklabels({'Int','DBX-','DBX+'});
xtickangle(45);
axis square;
box off;

% %%
%
% CorrVals = corr(FiringDBX);
% CorrVals = CorrVals - eye(size(CorrVals));
%
%
% subplot(4,4,13);
% molOrder = [find(CellOrder ==1)'; find(CellOrder ==2)'; find(CellOrder ==3)'];
% imagesc(CorrVals(molOrder,molOrder));
% colormap(colorcet('L1'));
% % set(gca,'XTick',1:size(order,1),'XTickLabel',{CellOrder(order)},'YTick',1:size(order,1),'YTickLabel',{CellOrder(order)});
% colorbar;
% axis square;
%
%
% [idx,c] = kmeans(CorrVals',2,'Distance','correlation','MaxIter',100);
% corrIds = idx;
% [~,order]=sort(idx);
%
% leadColor = 'r';
% lagColor = 'b';
%
%
% if trapz(mean(FiringDBX(1:50,corrIds==1),2)) > trapz(mean(FiringDBX(1:50,corrIds==2),2))
%     leadCorrIDs = find(corrIds==1);
%     lagCorrIDs = find(corrIds==2);
% else
%     leadCorrIDs = find(corrIds==2);
%     lagCorrIDs = find(corrIds==1);
% end
%
% subplot(4,4,14);
% FiringOrderedDBX_Trapz = [FiringDBX(:,leadCorrIDs),FiringDBX(:,lagCorrIDs)];
% imagesc(FiringOrderedDBX_Trapz');
% colormap(colorcet('L16'));
% % set(gca,'XTick',1:size(order,1),'XTickLabel',{CellOrder(order)},'YTick',1:size(order,1),'YTickLabel',{CellOrder(order)});
% colorbar;
% axis square;
% % title('Sorted by coorelations');
%
% subplot(4,4,15);
% shadedErrorBar(t,mean(FiringDBX(:,leadCorrIDs),2),std(FiringDBX(:,leadCorrIDs),[],2),'lineprops',{leadColor});
% hold on;
% shadedErrorBar(t,mean(FiringDBX(:,lagCorrIDs),2),std(FiringDBX(:,lagCorrIDs),[],2),'lineprops',{lagColor});
% xlabel('Time (sec)');
% ylabel('df/f');
% box off;
% set(gca, 'XLim',[-2,7]);
% axis square;
% title('Clusters from coorelation values');
%
% subplot(4,4,16)
% histogram(CellOrder(leadCorrIDs),'faceColor',leadColor)
% hold on ;
% histogram(CellOrder(lagCorrIDs),'faceColor',lagColor)
% box off;
% axis square;
% xticks([1,2,3]);
% xticklabels({'Int','vglut-','vglut+'})
% xtickangle(45);
% title('Type dist');


