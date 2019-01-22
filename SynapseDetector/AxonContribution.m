clear;

colors = cbrewer('qual','Paired',10);
startup

if ismac
    addpath(genpath('/Users/ashwin/Documents/'));
    df = readtable('/Users/ashwin/Documents/SynapseDetector/11252018.csv');
else
    addpath(genpath('/usr/people/ashwinv/seungmount/research/Ashwin/Scripts'));
    df = readtable('/usr/people/ashwinv/seungmount/research/Ashwin/SynapseDetector/11252018.csv');
end

%load('FiringRates.mat');
load('TAU.mat');
load('STA.mat');
load('AllCells.mat');
load('MelanieDBXCells.mat');


DbxCells = [76182,76183,76185,76186,76188,76189,76191,76199,76200];
OriginalCellOrderDBX = [8,9,10,11,21,12,15,2,3];
DbxTimeConstants = TAU(OriginalCellOrderDBX);

t = [-2:0.05:7]; % time period in sec

Firing = [STAall(:,OriginalCellOrderDBX) , DBX_vglut_neg/100 , DBX_vglut/100]; % divide by 100 to convert back from %
DBXpop = [DBX_vglut_neg/100 , DBX_vglut/100]; 
CellOrder = [ones(1,size(DbxCells,2)), 2* ones(1,size(DBX_vglut_neg,2)), 3* ones(1,size(DBX_vglut,2))];

% consider only the first 5 components of the trace
[A,B,C,D,E,F] = pca(Firing);
sprintf('first %d components capture %f of the data',5, sum(E(1:5)))
Firing = B(:,1:5)*A(:,1:5)'+ F;
normFiring = Firing(:,1:9);
%%  plot population data according to groups
% 2 --> Vglut_negative; blue
% 3 --> Vglut_positive; red

figure(1)
subplot(4,4,1)
shadedErrorBar(t,mean(Firing(:,CellOrder==2),2),std(Firing(:,CellOrder==2),[],2),'lineprops',{'b'});
hold on ; 
shadedErrorBar(t,mean(Firing(:,CellOrder==3),2),std(Firing(:,CellOrder==3),[],2),'lineprops',{'r'});
axis square;
title('Clusters from type');

%% Coorelation analyis of the cells.

CorrVals = corr(Firing);
[idx,c] = kmeans(Firing',2,'Distance','correlation');
corrIds = idx;
[~,order]=sort(idx);

leadColor = 'r';
lagColor = 'b';

figure(1);

if trapz(mean(Firing(1:50,corrIds==1),2)) > trapz(mean(Firing(1:50,corrIds==2),2))
    leadCorrIDs = find(corrIds==1);
    lagCorrIDs = find(corrIds==2);
else
    leadCorrIDs = find(corrIds==2);
    lagCorrIDs = find(corrIds==1);
end

subplot(4,4,5);
imagesc(CorrVals(order,order));
% set(gca,'XTick',1:size(order,1),'XTickLabel',{CellOrder(order)},'YTick',1:size(order,1),'YTickLabel',{CellOrder(order)});
colorbar;
axis square;
% title('Sorted by coorelations');

subplot(4,4,6);
shadedErrorBar(t,mean(Firing(:,leadCorrIDs),2),std(Firing(:,leadCorrIDs),[],2),'lineprops',{leadColor});
hold on;
plot(t,Firing(:,find(corrIds==1)<9),'Color',leadColor, 'LineWidth',0.5,'LineStyle','--');
shadedErrorBar(t,mean(Firing(:,lagCorrIDs),2),std(Firing(:,lagCorrIDs),[],2),'lineprops',{lagColor});
plot(t,Firing(:,find(corrIds==2)<9), 'Color',lagColor, 'LineWidth',0.5,'LineStyle','--');

box off;
set(gca, 'XLim',[-2,7]);
axis square;
title('Clusters from coorelation values');

leadNeurons = DbxCells(leadCorrIDs(leadCorrIDs<10));
lagNeurons  = DbxCells(lagCorrIDs(lagCorrIDs<10));

% find our cells in the population
subplot(4,4,7)
histogram(CellOrder(leadCorrIDs),'faceColor',leadColor)
hold on ;
histogram(CellOrder(lagCorrIDs),'faceColor',lagColor)
box off;
axis square;
title('Type dist');

% randomized cells in each group

subplot(4,4,8)
randomPop1 = randperm(160,80);
randomPop2 = 161-randomPop1;
shadedErrorBar(t,mean(Firing(:,randomPop1),2),std(Firing(:,randomPop1),[],2),'lineprops',{leadColor});
hold on;
shadedErrorBar(t,mean(Firing(:,randomPop2),2),std(Firing(:,randomPop2),[],2),'lineprops',{lagColor});
axis square;
title('Random clusters');

%% Extract Lead and Lag neuron partners from coorelation analysis.

% lead Neurons
Lead.PreSynapses = [];
Lead.PostSynapse = [];
Lead.PSDIDs = [];
Lead.PSDsize = [];
Lead.SynapseCount = [];
Lead.SynapseLocation = [];
Lead.SynapseLocationTransfromed = [];

for i = 1:numel(leadNeurons)
    [A,B] = SynapticPartners(leadNeurons(i),1,df);
    Lead(i).PreSynapses = A;
    Lead(i).PostSynapse = leadNeurons(i);
    Lead(i).PSDIDs = B;
    Lead(i).PSDsize = df.size(B);
    Lead(i).SynapseCount = histcounts(histcounts(Lead(i).PreSynapses, unique(Lead(i).PreSynapses)));
    Lead(i).SynapseLocation = PrePartnerCoordinates(B,df);
    Lead(i).SynapseLocationTransfromed = TransformPoints(Lead(i).SynapseLocation);
    clear A;
    clear B;
end

% lag neurons
Lag.PreSynapses = [];
Lag.PostSynapse = [];
Lag.PSDIDs = [];
Lag.PSDsize = [];
Lag.SynapseCount = [];
Lag.SynapseLocation = [];
Lag.SynapseLocationTransfromed = [];

for i = 1:numel(lagNeurons)
    [A,B] = SynapticPartners(lagNeurons(i),1,df);
    Lag(i).PreSynapses = A;
    Lag(i).PostSynapse = lagNeurons(i);
    Lag(i).PSDIDs = B;
    Lag(i).PSDsize = df.size(B);
    Lag(i).SynapseCount = histcounts(histcounts(Lag(i).PreSynapses,unique(Lag(i).PreSynapses)));
    Lag(i).SynapseLocation = PrePartnerCoordinates(B,df);
    Lag(i).SynapseLocationTransfromed = TransformPoints(Lag(i).SynapseLocation);
    clear A;
    clear B;
end

figure(2)

subplot(4,4,1)

histogram(vertcat(Lead.PSDsize)/sum(vertcat(Lead.PSDsize)),'FaceColor',leadColor,'EdgeColor','none');
hold on;
histogram(vertcat(Lag.PSDsize)/sum(vertcat(Lag.PSDsize)),'FaceColor',lagColor,'EdgeColor','none');
box off;
axis square;


leadAxons = vertcat(Lead.PreSynapses);
numberOfLeadAxons = numel(leadAxons);
leadAxons = leadAxons(leadAxons<1e5);

lagAxons = vertcat(Lag.PreSynapses);
numberOfLagAxons = numel(lagAxons);
lagAxons = lagAxons(lagAxons<1e5);

subplot(4,4,2)
histogram(histcounts(vertcat(Lead.PreSynapses),unique(vertcat(Lead.PreSynapses))),'FaceColor',leadColor,'EdgeColor','none');
hold on;
histogram(histcounts(vertcat(Lag.PreSynapses),unique(vertcat(Lag.PreSynapses))),'FaceColor',lagColor,'EdgeColor','none');
box off;
axis square;


%% plot the axons on the Z-brain atlas
figure(3);

subplot(1,3,1)
transform_swc_AV(unique(leadAxons)',colors(5,:),colors(6,:),[],false);
subplot(1,3,2);
transform_swc_AV(unique(lagAxons)',colors(1,:),colors(2,:),[],false);


% plot the location of the synapses

leadSynapseLocations = vertcat(Lead.SynapseLocationTransfromed);
lagSynapseLocations = vertcat(Lag.SynapseLocationTransfromed);

subplot(1,3,3);
transform_swc_AV(leadNeurons,colors(5,:),colors(6,:),[],false);
hold on;
transform_swc_AV(lagNeurons,colors(1,:),colors(2,:),[],false);
scatter3(leadSynapseLocations(:,1), leadSynapseLocations(:,2),leadSynapseLocations(:,3),...
    8,'o','MarkerFaceColor',leadColor,'MarkerEdgeColor','none');
hold on;
scatter3(lagSynapseLocations(:,1), lagSynapseLocations(:,2),lagSynapseLocations(:,3),...
    8,'o','MarkerFaceColor',lagColor,'MarkerEdgeColor','none');
daspect([1,1,1]);
set(gca, 'BoxStyle','full','YDir','reverse','ZDir','reverse');
set(gca,'ZLim',[0 276],'XLim',[0 495.558], 'YLim',[0 1121.988]);
%% Extract lead/lag axons identity

for i =1:numel(leadNeurons)
    for j = 1:numel(Lead(i).PreSynapses)
        l = isSaccade(Lead(i).PreSynapses(j));
        Lead(i).isSaccade(j,1) = l;
        if l == 1
            Lead(i).saccadicSynapses(j) = df.psd_segid(df.presyn_segid == Lead(i).PreSynapses(j) ...
                & df.postsyn_segid == leadNeurons(i));
        else 
            continue
        end
    end
end

leadSaccade = leadSaccade(find(leadSaccade==1),2);

for i =1:numel(lagNeurons)
    for j = 1:numel(Lag(i).PreSynapses)
        l = isSaccade(Lag(i).PreSynapses(j));
        Lag(i).isSaccade(j,1) = l;
    end
end
lagSaccade = lagSaccade(find(lagSaccade==1),2);




%%

corrAndSizeDBX = [];
for i = 1:length(DbxCells)
    for j = 1:length(DbxCells)
        if i == j
            corrAndSizeDBX = [corrAndSizeDBX];
        else
            [preSynapticInputs1,preSynapticInputs1PSD] = SynapticPartners(DbxCells(i),1,df);
            [preSynapticInputs2,preSynapticInputs2PSD] = SynapticPartners(DbxCells(j),1,df);
            [commonInputs,loc1,loc2] = intersect(preSynapticInputs1, preSynapticInputs2);
            psdSize1 = df.size(preSynapticInputs1PSD(loc1));
            psdSize2 = df.size(preSynapticInputs2PSD(loc2));
            tempCorr = corrcoef(normFiring(:,OriginalCellOrderDBX(i)), normFiring(:,OriginalCellOrderDBX(j)));
            corrAndSizeDBX = [corrAndSizeDBX; commonInputs, tempCorr(1,2).*ones(size(psdSize1,1),1), psdSize1, psdSize2];
            
            clear {'preSynapticInputs1','preSynapticInputs1PSD', 'psdSize1','loc1'};
            clear {'preSynapticInputs2', 'preSynapticInputs2PSD','psdSize2', 'loc2'};
            clear commonInputs;
        end
    end
end


uniqueCommonInputs = unique(corrAndSizeDBX(:,1));
index  = 1;

for i =1:size(uniqueCommonInputs,1)
    [A,B] = SynapticPartners(uniqueCommonInputs(i),2,df);
    [lia,lib] = ismember(DbxCells,A);
    if sum(lia)>1
        averagedFiring(index,:) = mean(normFiring(:,OriginalCellOrderDBX(find(lia))),2);
        axon(index).ID  = uniqueCommonInputs(i);
        index = index+1;
    end
end

clear lia;
clear lib;

% averageFiringDBX = mean(zeroedPop,2);
% 
% for i = 1:size(axon,2)
%    % subplot(12,12,i)
%    % plot(t,(averagedFiring(i,:)'- averageFiringDBX));
%     DiffFromAverage(i,:) = (averagedFiring(i,:)'- averageFiringDBX);
%     line([-2,7], [0,0],'color','k');
%     XLim = [-2,7];
%     box off;
%     axis square;
%     axis off;
%     % title(axonID(i));
% end


for i = 1:size(axon,2)
    [postCells,postCellsPSD] = SynapticPartners(axon(i).ID,2,df);
    [axon(i).PostDBX, lia, lib] = intersect(postCells,DbxCells);
    axon(i).PostPSD = postCellsPSD(lia);
    for j = 1:size(lia,1)
        axon(i).PostPSD_x(j) = df.centroid_x(df.psd_segid == axon(i).PostPSD(j));
        axon(i).PostPSD_y(j) = df.centroid_y(df.psd_segid == axon(i).PostPSD(j));
        axon(i).PostPSD_z(j) = df.centroid_z(df.psd_segid == axon(i).PostPSD(j));
    end
end


%% sort the traces
%clear DiffFromAverage;

DiffFromAverage = averagedFiring;

clear clust
for i = 1:1:30
    clust(:,i) = kmeans(DiffFromAverage,i,'Distance','correlation');
end

eva = evalclusters(DiffFromAverage,clust,'silhouette');
clear idx;
clear c;

[idx,c] = kmeans(DiffFromAverage,2,'Distance','correlation','MaxIter',100);


if size(find(idx == 1),1) > size(find(idx == 2))
    lead = 'r';
    lag = 'b';
else
    lead = 'b'
    lag = 'r'
end

cluster1 = find(idx==1);
subplot(4,4,5)
shadedErrorBar(t, mean(DiffFromAverage(cluster1,:)),std(DiffFromAverage(cluster1,:)),{lead},0.2);
%legend({'culster1'})
hold on;
coordX = [];
coordY = [];
coordZ = [];

for i = 1:size(cluster1,1)
    coordX = [coordX, axon(cluster1(i)).PostPSD_x];
    coordY = [coordY, axon(cluster1(i)).PostPSD_y];
    coordZ = [coordZ, axon(cluster1(i)).PostPSD_z];
end


subplot(4,4,6)
cluster1Coords = [coordX', coordY', coordZ'];
cluster1CoordsTransformed = TransformPoints(cluster1Coords);

scatter3(cluster1CoordsTransformed(:,1),cluster1CoordsTransformed(:,2), cluster1CoordsTransformed(:,3), 'MarkerFaceColor',lead,'MarkerFaceAlpha',0.2);
hold on;
scatter3(mean(cluster1CoordsTransformed(:,1)), mean(cluster1CoordsTransformed(:,2)), mean(cluster1CoordsTransformed(:,3)), 100, 'MarkerFaceColor',lead);


cluster2 = find(idx == 2);

subplot(4,4,5)
shadedErrorBar(t, mean(DiffFromAverage(cluster2,:)),std(DiffFromAverage(cluster2,:)),{lag},0.2);
box off;
set(gca,'XLim',[-2,7]);
title('Observed mean activity profiles');

coordX = [];
coordY = [];
coordZ = []; 

for i = 1:size(cluster2,1)
    coordX = [coordX, axon(cluster2(i)).PostPSD_x];
    coordY = [coordY, axon(cluster2(i)).PostPSD_y];
    coordZ = [coordZ, axon(cluster2(i)).PostPSD_z];
end

cluster2Coords = [coordX', coordY', coordZ'];
cluster2CoordsTransformed = TransformPoints(cluster2Coords);

subplot(4,4,6)
scatter3(cluster2CoordsTransformed(:,1),cluster2CoordsTransformed(:,2), cluster2CoordsTransformed(:,3), 'MarkerFaceColor',lag, 'MarkerFaceAlpha',0.2);
scatter3(mean(cluster2CoordsTransformed(:,1)), mean(cluster2CoordsTransformed(:,2)), mean(cluster2CoordsTransformed(:,3)), 100, 'MarkerFaceColor',lag);
title('Synapse positions');
%% boot strap with the population data

for i = 1:10000
    randomizeCellOrder = randperm(size(DBXpop,2),randperm(6,1)); % randomly select 1-5 neurons from the population
    randomSampling(i,:) = mean(DBXpop(:,randomizeCellOrder),2);
   % randomSampleingAverageDiff(i,:) = randomSampling(i,:) - averageFiringDBX' ; 
end

% for i = 1:1:10
%     clust(:,i) = kmeans(randomSampleingAverageDiff,i+1,'Distance','cosine');
% end
% 
% eva = evalclusters(randomSampleingAverageDiff,clust,'silhouette');
clear idx;
clear c;
[idx,c] = kmeans(randomSampling,2,'Distance','cosine','MaxIter',100);

if size(find(idx == 1),1) > size(find(idx == 2))
    lead = 'r';
    lag = 'b';
else
    lead = 'b'
    lag = 'r'
end


subplot(4,4,9)
randomcluster1 = find(idx == 1);
shadedErrorBar(t, mean(randomSampling(randomcluster1,:)),std(randomSampling(randomcluster1,:)),{lead,'LineStyle','-'},0.2);
hold on;
randomcluster2 = find(idx == 2);
shadedErrorBar(t, mean(randomSampling(randomcluster2,:)),std(randomSampling(randomcluster2,:)),{lag,'LineStyle','-'},0.2)
box off;
set(gca, 'XLim',[-2,7]);
title('Randomized mean activity profiles');

%% post synaptic partners in each cluster


cluster1PostIDs = [];
for i =1:size(cluster1,1)
    cluster1PostIDs = [cluster1PostIDs;axon(cluster1(i)).PostDBX];
end

cluster2PostIDs = [];
for i =1:size(cluster2,1)
    cluster2PostIDs = [cluster2PostIDs;axon(cluster2(i)).PostDBX];
end

subplot(4,4,13)
 hist1  = histcounts(cluster1PostIDs,'NumBins',size(unique(cluster1PostIDs),1));
 %hist1  = histogram(cluster1PostIDs,'NumBins',size(unique(cluster1PostIDs),1),'FaceColor','r','FaceAlpha',0.2);
% hold on;
 hist2 = histcounts(cluster2PostIDs,'NumBins',size(unique(cluster2PostIDs),1));
 %hist2 = histogram(cluster2PostIDs,'NumBins',size(unique(cluster2PostIDs),1),'FaceColor','b','FaceAlpha',0.2);

% box off
% title('Most represented neuron');

% subplot(4,4,14)
[~,loc] = sort(hist1);
[clust1max,~] = SynapticPartners(DbxCells(loc(end)),1,df);
[clust1penmax,~] = SynapticPartners(DbxCells(loc(end-1)),1,df);
[cluster1mostDominiant,~] = intersect([clust1max;clust1penmax],uniqueCommonInputs);

plot(t,Firing(:,OriginalCellOrderDBX(loc(end))),lead);
hold on;
plot(t,Firing(:,OriginalCellOrderDBX(loc(end-1))),lead);

[~,loc] = sort(hist2);
[clust2max,~] = SynapticPartners(DbxCells(loc(end)),1,df);
[clust2penmax,~] = SynapticPartners(DbxCells(loc(end-1)),1,df);
[cluster2mostDominiant,~] = intersect([clust2max;clust2penmax],uniqueCommonInputs);

plot(t,Firing(:,OriginalCellOrderDBX(loc(end))),lag);
hold on;
plot(t,Firing(:,OriginalCellOrderDBX(loc(end-1))),lag);
box off;

%% plot the clusters on the Zbrian atlas

subplot(4,4,14)
transform_swc_AV(cluster1mostDominiant',colors(1,:),colors(2,:),[],false);
subplot(4,4,15);
transform_swc_AV(cluster2mostDominiant',colors(3,:),colors(4,:),[],false);

subplot(4,4,16)
transform_swc_AV(setdiff(cluster1mostDominiant,cluster2mostDominiant)',colors(5,:),colors(6,:),[],false);

%%
for i = 1:9
    [DBXPre, DBXPrePSDid] = SynapticPartners(DbxCells(i),1,df);
     DBXPreCell(i).ID = DbxCells(i);
     [DBXPreCell(i).Cluster1,lib1] = intersect([axon(cluster1).ID]',DBXPre);
     [DBXPreCell(i).Cluster2,lib2] = intersect([axon(cluster2).ID]',DBXPre);
     DBXPreCell(i).X =[];
     DBXPreCell(i).Y =[];
     DBXPreCell(i).Z =[];
     for j = 1:size(lib1,1)
        DBXPreCell(i).X = [DBXPreCell(i).X; df.centroid_x(df.psd_segid == DBXPrePSDid(lib1(i)))];
        DBXPreCell(i).Y = [DBXPreCell(i).Y; df.centroid_y(df.psd_segid == DBXPrePSDid(lib1(i)))];
        DBXPreCell(i).Z = [DBXPreCell(i).Z; df.centroid_z(df.psd_segid == DBXPrePSDid(lib1(i)))];
     end   
    % DBXCluster1{i,2} = [axon(DBXCluster1{i,1}).PostPSD_x, axon(DBXCluster1{i,1}).PostPSD_y, axon(DBXCluster1{i,1}).PostPSD_z];
    % DBXCluster2{i,2} = [axon(DBXCluster2{i,1}).PostPSD_x, axon(DBXCluster2{i,1}).PostPSD_y, axon(DBXCluster2{i,1}).PostPSD_z];
end

%% shuffle the donwstream partner of the axon.

randomDifferenceFromAverage = zeros(140,181,500);
for j = 1:100
   % OriginalCellOrderDBXPerumted = OriginalCellOrderDBX(randperm(length(OriginalCellOrderDBX)));
    OriginalCellOrderDBXPerumted = randperm(size(DBXpop,2),9);
    normFiringPermuted = DBXpop(:,OriginalCellOrderDBXPerumted);
    index  = 1;
    for i =1:size(uniqueCommonInputs,1)
        [A,B] = SynapticPartners(uniqueCommonInputs(i),2,df);
        [lia,lib] = ismember(DbxCells,A);
        if sum(lia)>1
            averagedFiring(index,:) = mean(normFiringPermuted(:,find(lia)),2);
            randomDifferenceFromAverage(index,:,j) = averagedFiring(index,:) - averageFiringDBX(:,1)';
            axon(index).ID  = uniqueCommonInputs(i);
            index = index+1;
        end
    end
    disp(j);
end
%%
for i = 1:100
    randomIndex = randperm(100,1);
    temp  = randomDifferenceFromAverage(:,:,randomIndex);
    [idx,c] = kmeans(temp,2,'Distance','cosine','MaxIter',100);
    randomcluster1 = find(idx == 1);
    temp1 = mean(temp(randomcluster1,:));
    randomcluster2 = find(idx == 2);
    temp2 = mean(temp(randomcluster2,:));
    clear temp;
    if trapz(temp1(1:50)) > 0
        meanCluster1(i,:) = temp2;
        meanCluster2(i,:) = temp1;
    else
        meanCluster1(i,:) = temp1;
        meanCluster2(i,:) = temp2;
    end
end


figure(2);
shadedErrorBar(t, mean(meanCluster1),std(meanCluster1),{'g','lineStyle','--'},0.2);
hold on;
shadedErrorBar(t, mean(meanCluster2),std(meanCluster2),{'c', 'lineStyle','--'},0.2)
