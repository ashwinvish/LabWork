clear;
addpath(genpath('/Users/ashwin/Documents/'));

colors = cbrewer('qual','Dark2',10);
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
load('MelanieALXCells.mat');



DbxCells = [76182,76183,76185,76186,76188,76189,76191,76199,76200];
OriginalCellOrderDBX = [8,9,10,11,21,12,15,2,3];
DbxTimeConstants = TAU(OriginalCellOrderDBX);

t = [-2:0.05:7]; % time period in sec

Firing = [STAall(:,OriginalCellOrderDBX) , DBX_vglut_neg/100 , DBX_vglut/100]; % divide by 100 to convert back from %
DBXpop = [DBX_vglut_neg/100 , DBX_vglut/100]; % melanies data is reported as %
CellOrder = [ones(1,size(DbxCells,2)), 2* ones(1,size(DBX_vglut_neg,2)), 3* ones(1,size(DBX_vglut,2))];

% consider only the first 5 components of the trace
[A,B,C,D,E,F] = pca(Firing);
sprintf('first %d components capture %f of the data',5, sum(E(1:5)))
Firing = B(:,1:5)*A(:,1:5)'+ F;
normFiring = Firing(:,1:9);


temp1 = colorcet('CBTL1','N',5); %  linear-tritanopic_krjcw_5-98_c46_n256
temp2 = colorcet('CBL2','N',5);  %  linear-protanopic-deuteranopic_kbw_5-98_c40_n256

leadColor = temp1(3,:); % dark red, CB safe
lagColor = temp2(3,:); % dark blue, CB safe

clear temp1;
clear temp2;

PartnerColors = colorcet('CBTD1','N',10);
lightRed = PartnerColors(10,:);
lightBlue = PartnerColors(1,:);


%%  plot population data according to groups
% 2 --> Vglut_negative; blue
% 3 --> Vglut_positive; red

figure(1)
subplot(4,4,1)
shadedErrorBar(t,mean(Firing(:,CellOrder==2),2),std(Firing(:,CellOrder==2),[],2),'lineprops',{'Color',lagColor});
hold on ;
shadedErrorBar(t,mean(Firing(:,CellOrder==3),2),std(Firing(:,CellOrder==3),[],2),'lineprops',{'Color',leadColor});
xlabel('Time (sec)');
ylabel('df/f');
axis square;

title('Saccade triggered fluorescence (by labels)');

molOrder = [find(CellOrder ==1)'; find(CellOrder ==2)'; find(CellOrder ==3)'];

subplot(4,4,3)
imagesc(Firing(:,molOrder)');
colormap(colorcet('L16'));
colorbar;
set(gca,'XTick',[0,41,81,121,161],'XTickLabel',[-2,0,2,4,6]);
xlabel('Time(sec)');
box off;
ylabel('Neurons');
daspect([1,1,1]);

%% Coorelation analyis of the cells.
% CorrVals = corr(Firing);
% CorrVals = CorrVals - eye(size(CorrVals));
% cutoff = 0.5;
% Z = linkage(CorrVals,'complete','correlation')

eva = evalclusters(Firing','kmeans','silhouette','Distance','correlation','KList',1:6);

subplot(4,4,4)
plot(eva);
title('Optimal clusters');
axis square;
box off;

% subplot(4,4,2);
% molOrder = [find(CellOrder ==1)'; find(CellOrder ==2)'; find(CellOrder ==3)'];
% imagesc(CorrVals(molOrder,molOrder));
% colormap(colorcet('L1'));
% % set(gca,'XTick',1:size(order,1),'XTickLabel',{CellOrder(order)},'YTick',1:size(order,1),'YTickLabel',{CellOrder(order)});
% colorbar;
% axis square;


% [idx,c] = kmeans(CorrVals',2,'Distance','correlation','MaxIter',100);
% corrIds = idx;
% [~,order]=sort(idx);


figure(1);

if max(mean(Firing(1:50,eva.OptimalY==1),2)) > max(mean(Firing(1:50,eva.OptimalY==2),2))
    leadCorrIDs = find(eva.OptimalY==1);
    lagCorrIDs = find(eva.OptimalY==2);
else
    leadCorrIDs = find(eva.OptimalY==2);
    lagCorrIDs = find(eva.OptimalY==1);
end

order = [find(eva.OptimalY==1);find(eva.OptimalY ==2)];

subplot(4,4,5);
imagesc(Firing(:,order)');
colormap(colorcet('L16'));
colorbar;
set(gca,'XTick',[0,41,81,121,161],'XTickLabel',[-2,0,2,4,6]);
xlabel('Time(sec)');
box off;
ylabel('Neurons');
daspect([1,1,1]);


subplot(4,4,6);
shadedErrorBar(t,mean(Firing(:,leadCorrIDs),2),std(Firing(:,leadCorrIDs),[],2),'lineProps',{'Color',leadColor});
hold on;
shadedErrorBar(t,mean(Firing(:,lagCorrIDs),2),std(Firing(:,lagCorrIDs),[],2),'lineprops',{'Color',lagColor});
xlabel('Time (sec)');
ylabel('df/f');
box off;
set(gca, 'XLim',[-2,7]);
axis square;
title('Clusters from coorelation values');

subplot(4,4,7)
plot(t,Firing(:,leadCorrIDs(find(leadCorrIDs<=9))),'Color',leadColor, 'LineWidth',1);
xlabel('Time (sec)');
ylabel('df/f');
box off;
set(gca, 'XLim',[-2,7]);
axis square;
title('Clusters from coorelation values');


subplot(4,4,8)
plot(t,Firing(:,lagCorrIDs(find(lagCorrIDs<=9))), 'Color',lagColor, 'LineWidth',1);
xlabel('Time (sec)');
ylabel('df/f');
box off;
set(gca, 'XLim',[-2,7]);
axis square;
title('Clusters from coorelation values');


leadNeurons = DbxCells(leadCorrIDs(leadCorrIDs<10));
lagNeurons  = DbxCells(lagCorrIDs(lagCorrIDs<10));

save('leadNeurons.mat','leadNeurons');
save('lagNeurons.mat','lagNeurons');

% find our cells in the population
subplot(4,4,9)
histogram(CellOrder(leadCorrIDs),'faceColor',leadColor)
hold on ;
histogram(CellOrder(lagCorrIDs),'faceColor',lagColor)
box off;
axis square;
xticks([1,2,3]);
xticklabels({'Int','vglut-','vglut+'})
xtickangle(45);
title('Type dist');

% randomized cells in each group

subplot(4,4,10)
randomPop1 = randperm(160,80);
randomPop2 = 161-randomPop1;
shadedErrorBar(t,mean(Firing(:,randomPop1),2),std(Firing(:,randomPop1),[],2),'lineprops',{'Color',leadColor});
hold on;
shadedErrorBar(t,mean(Firing(:,randomPop2),2),std(Firing(:,randomPop2),[],2),'lineprops',{'Color',lagColor});
axis square;
xlabel('Time (sec)');
ylabel('df/f');
title('Random clusters');

%%

load smallABDneurons.mat
load largeABDneurons.mat

for i = 1:numel(DbxCells)
    DBX(i) = InputsByClass(DbxCells(i),df);
end

%% Saccadic module

UniqueSaccadicAxons = unique(vertcat(DBX.Saccadic));
    ix1 = 1;
    ix2 = 1;
    
for i = 1:numel(UniqueSaccadicAxons)
    [A,~] = SynapticPartners(UniqueSaccadicAxons(i),2,df);
    A = A(A<1e5);
    leadSum  = sum(ismember(A,leadNeurons));
    lagSum = sum(ismember(A,lagNeurons));
    smallCount = sum(ismember(A,smallABDCellIDs));
    largeCount = sum(ismember(A,largeABDCellIDs));

    if leadSum>lagSum
        Lead(ix1).SaccadeAxonID = UniqueSaccadicAxons(i);
        Lead(ix1).SaccadeDiff = leadSum-lagSum;
        Lead(ix1).SaccadeMotorDist = MotorDiff(UniqueSaccadicAxons(i),df);
        Lead(ix1).SaccadeSmallCount = smallCount;
        Lead(ix1).SaccadeLargeCount = largeCount;
        Lead(ix1).SaccadeLargeDiffSmall = largeCount-smallCount;
        ix1 = ix1+1;
    else
        Lag(ix2).SaccadeAxonID = UniqueSaccadicAxons(i);
        Lag(ix2).SaccadeDiff = leadSum-lagSum;
        Lag(ix2).SaccadeMotorDist = MotorDiff(UniqueSaccadicAxons(i),df);
        Lag(ix2).SaccadeSmallCount = smallCount;
        Lag(ix2).SaccadeLargeCount = largeCount;
        Lag(ix2).SaccadeLargeDiffSmall = largeCount-smallCount;
        ix2 = ix2+1;
    end
    clear A;
end

for i = 1:numel(DbxCells)
    tempLeadIDs = find(ismember(DBX(i).Inputs, [Lead.SaccadeAxonID]));
    tempLagIDs = find(ismember(DBX(i).Inputs, [Lag.SaccadeAxonID]));
    DBX(i).LeadSaccadicPathLength = DBX(i).PathLength(tempLeadIDs)/max(Pvec_tree(DBX(i).Tree{1}));
    DBX(i).LeadSaccadicPSD = DBX(i).PSDsize(tempLeadIDs);
    LeadSaccadicPathLengths(i,:) = histcounts(DBX(i).LeadSaccadicPathLength,0:0.1:1);
    
    DBX(i).LagSaccadicPathLength = DBX(i).PathLength(tempLagIDs)/max(Pvec_tree(DBX(i).Tree{1}));
    DBX(i).LagSaccadicPSD = DBX(i).PSDsize(tempLagIDs);
    LagSaccadicPathLengths(i,:) = histcounts(DBX(i).LagSaccadicPathLength,0:0.1:1);    
    
    clear tempLeadIDs;
    clear tempLagIDs;
end

subplot(4,4,1)
scatter([Lead.SaccadeDiff],[Lead.SaccadeMotorDist],20,'filled','MarkerFaceColor',leadColor);
hold on;
scatter([Lag.SaccadeDiff],[Lag.SaccadeMotorDist],20,'filled','MarkerFaceColor',lagColor);
line([-10,10],[0,0],'color','k','LineStyle','--');
line([0,0],[-100,100],'color','k','LineStyle','--');
axis square;
xlabel('Synapse diff');
ylabel('A-Ai diff');

subplot(4,4,2)
histogram(vertcat(DBX.LeadSaccadicPathLength),10,'FaceColor',leadColor,'EdgeColor','none');
hold on;
histogram(vertcat(DBX.LagSaccadicPathLength),10,'FaceColor',lagColor,'EdgeColor','none');
axis square;
box off;
xlabel('Norm. pathlength');
ylabel('count')

subplot(4,4,3)
shadedErrorBar(0.1:0.1:1,mean(LeadSaccadicPathLengths),std(LeadSaccadicPathLengths)/sqrt(9),'lineProps',{'Color',leadColor,'LineWidth',2});
hold on;
shadedErrorBar(0.1:0.1:1,mean(LagSaccadicPathLengths),std(LagSaccadicPathLengths)/sqrt(9),'lineProps',{'Color',lagColor,'LineWidth',2});
axis square;
xlabel('Norm pathlength');
ylabel('count per cell');

%% Vestibular Module


UniqueVestibularAxons = unique(vertcat(DBX.Vestibular));
    ix1 = 1;
    ix2 = 1;  
for i = 1:numel(UniqueVestibularAxons)
    [A,~] = SynapticPartners(UniqueVestibularAxons(i),2,df);
    %A = A(A<1e5);
    leadSum  = sum(ismember(A,leadNeurons));
    lagSum = sum(ismember(A,lagNeurons));
    smallCount = sum(ismember(A,smallABDCellIDs));
    largeCount = sum(ismember(A,largeABDCellIDs));

    if leadSum>lagSum
        Lead(ix1).VestibularAxonID = UniqueVestibularAxons(i);
        Lead(ix1).VestibularDiff = leadSum-lagSum;
        Lead(ix1).VestibularMotorDist = MotorDiff(UniqueVestibularAxons(i),df);
        Lead(ix1).VestibularSmallCount = smallCount;
        Lead(ix1).VestibularLargeCount = largeCount;
        Lead(ix1).VestibularLargeDiffSmall = largeCount-smallCount;
        ix1 = ix1+1;
    else
        Lag(ix2).VestibularAxonID = UniqueVestibularAxons(i);
        Lag(ix2).VestibularDiff = leadSum-lagSum;
        Lag(ix2).VestibularMotorDist = MotorDiff(UniqueVestibularAxons(i),df);
        Lag(ix2).VestibularSmallCount = smallCount;
        Lag(ix2).VestibularLargeCount = largeCount;
        Lag(ix2).VestibularLargeDiffSmall = largeCount-smallCount;
        ix2 = ix2+1;
    end
    clear A;
end

for i = 1:numel(DbxCells)
    tempLeadIDs = find(ismember(DBX(i).Inputs, [Lead.VestibularAxonID]));
    tempLagIDs = find(ismember(DBX(i).Inputs, [Lag.VestibularAxonID]));
    DBX(i).LeadVestibularPathLength = DBX(i).PathLength(tempLeadIDs)/max(Pvec_tree(DBX(i).Tree{1}));
    DBX(i).LeadVestibularPSD = DBX(i).PSDsize(tempLeadIDs);
    LeadVestibularPathLengths(i,:) = histcounts(DBX(i).LeadVestibularPathLength,0:0.1:1);
    
    DBX(i).LagVestibularPathLength = DBX(i).PathLength(tempLagIDs)/max(Pvec_tree(DBX(i).Tree{1}));
    DBX(i).LagVestibularPSD = DBX(i).PSDsize(tempLagIDs);
    LagVestibularPathLengths(i,:) = histcounts(DBX(i).LagVestibularPathLength,0:0.1:1);    
    
    clear tempLeadIDs;
    clear tempLagIDs;
end


subplot(4,4,5)
scatter([Lead.VestibularDiff],[Lead.VestibularMotorDist],20,'filled','MarkerFaceColor',leadColor);
hold on;
scatter([Lag.VestibularDiff],[Lag.VestibularMotorDist],20,'filled','MarkerFaceColor',lagColor);
line([-10,10],[0,0],'color','k','LineStyle','--');
line([0,0],[-100,100],'color','k','LineStyle','--');
axis square;
xlabel('Synapse diff');
ylabel('A-Ai diff');

subplot(4,4,6)
histogram(vertcat(DBX.LeadVestibularPathLength),10,'FaceColor',leadColor,'EdgeColor','none');
hold on;
histogram(vertcat(DBX.LagVestibularPathLength),10,'FaceColor',lagColor,'EdgeColor','none');
axis square;
box off;
xlabel('Norm. pathlength');
ylabel('count')

subplot(4,4,7)
shadedErrorBar(0.1:0.1:1,mean(LeadVestibularPathLengths),std(LeadVestibularPathLengths)/sqrt(9),'lineProps',{'Color',leadColor,'LineWidth',2});
hold on;
shadedErrorBar(0.1:0.1:1,mean(LagVestibularPathLengths),std(LagVestibularPathLengths)/sqrt(9),'lineProps',{'Color',lagColor,'LineWidth',2});
axis square;
xlabel('Norm pathlength');
ylabel('count per cell');

%% Integrator Module


UniqueIntegratorAxons = unique(vertcat(DBX.Integrator));
    ix1 = 1;
    ix2 = 1;  
for i = 1:numel(UniqueIntegratorAxons)
    [A,~] = SynapticPartners(UniqueIntegratorAxons(i),2,df);
    %A = A(A<1e5);
    leadSum  = sum(ismember(A,leadNeurons));
    lagSum = sum(ismember(A,lagNeurons));
    smallCount = sum(ismember(A,smallABDCellIDs));
    largeCount = sum(ismember(A,largeABDCellIDs));

    if leadSum>lagSum
        Lead(ix1).IntegratorAxonID = UniqueIntegratorAxons(i);
        Lead(ix1).IntegratorDiff = leadSum-lagSum;
        Lead(ix1).IntegratorMotorDist = MotorDiff(UniqueIntegratorAxons(i),df);
        Lead(ix1).IntegratorSmallCount = smallCount;
        Lead(ix1).IntegratorLargeCount = largeCount;
        Lead(ix1).IntegratorLargeDiffSmall = largeCount-smallCount;
        ix1 = ix1+1;
    else
        Lag(ix2).IntegratorAxonID = UniqueIntegratorAxons(i);
        Lag(ix2).IntegratorDiff = leadSum-lagSum;
        Lag(ix2).IntegratorMotorDist = MotorDiff(UniqueIntegratorAxons(i),df);
        Lag(ix2).IntegratorSmallCount = smallCount;
        Lag(ix2).IntegratorLargeCount = largeCount;
        Lag(ix2).IntegratorLargeDiffSmall = largeCount-smallCount;
        ix2 = ix2+1;
    end
    clear A;
end

for i = 1:numel(DbxCells)
    tempLeadIDs = find(ismember(DBX(i).Inputs, [Lead.IntegratorAxonID]));
    tempLagIDs = find(ismember(DBX(i).Inputs, [Lag.IntegratorAxonID]));
    DBX(i).LeadIntegratorPathLength = DBX(i).PathLength(tempLeadIDs)/max(Pvec_tree(DBX(i).Tree{1}));
    DBX(i).LeadIntegratorPSD = DBX(i).PSDsize(tempLeadIDs);
    LeadIntegratorPathLengths(i,:) = histcounts(DBX(i).LeadIntegratorPathLength,0:0.1:1);
    
    DBX(i).LagIntegratorPathLength = DBX(i).PathLength(tempLagIDs)/max(Pvec_tree(DBX(i).Tree{1}));
    DBX(i).LagIntegratorPSD = DBX(i).PSDsize(tempLagIDs);
    LagIntegratorPathLengths(i,:) = histcounts(DBX(i).LagIntegratorPathLength,0:0.1:1);    
    
    clear tempLeadIDs;
    clear tempLagIDs;
end


subplot(4,4,9)
scatter([Lead.IntegratorDiff],[Lead.IntegratorMotorDist],20,'filled','MarkerFaceColor',leadColor);
hold on;
scatter([Lag.IntegratorDiff],[Lag.IntegratorMotorDist],20,'filled','MarkerFaceColor',lagColor);
line([-10,10],[0,0],'color','k','LineStyle','--');
line([0,0],[-100,100],'color','k','LineStyle','--');
axis square;
xlabel('Synapse diff');
ylabel('A-Ai diff');

subplot(4,4,10)
histogram(vertcat(DBX.LeadIntegratorPathLength),10,'FaceColor',leadColor,'EdgeColor','none');
hold on;
histogram(vertcat(DBX.LagIntegratorPathLength),10,'FaceColor',lagColor,'EdgeColor','none');
axis square;
box off;
xlabel('Norm. pathlength');
ylabel('count')

subplot(4,4,11)
shadedErrorBar(0.1:0.1:1,mean(LeadIntegratorPathLengths),std(LeadIntegratorPathLengths)/sqrt(9),'lineProps',{'Color',leadColor,'LineWidth',2});
hold on;
shadedErrorBar(0.1:0.1:1,mean(LagIntegratorPathLengths),std(LagIntegratorPathLengths)/sqrt(9),'lineProps',{'Color',lagColor,'LineWidth',2});
axis square;
xlabel('Norm pathlength');
ylabel('count per cell');

%% Contra Module

UniqueContraAxons = unique(vertcat(DBX.Contra));
    ix1 = 1;
    ix2 = 1;  
for i = 1:numel(UniqueContraAxons)
    [A,~] = SynapticPartners(UniqueContraAxons(i),2,df);
    %A = A(A<1e5);
    leadSum  = sum(ismember(A,leadNeurons));
    lagSum = sum(ismember(A,lagNeurons));
    smallCount = sum(ismember(A,smallABDCellIDs));
    largeCount = sum(ismember(A,largeABDCellIDs));

    if leadSum>lagSum
        Lead(ix1).ContraAxonID = UniqueContraAxons(i);
        Lead(ix1).ContraDiff = leadSum-lagSum;
        Lead(ix1).ContraMotorDist = MotorDiff(UniqueContraAxons(i),df);
        Lead(ix1).ContraSmallCount = smallCount;
        Lead(ix1).ContraLargeCount = largeCount;
        Lead(ix1).ContraLargeDiffSmall = largeCount-smallCount;
        ix1 = ix1+1;
    else
        Lag(ix2).ContraAxonID = UniqueContraAxons(i);
        Lag(ix2).ContraDiff = leadSum-lagSum;
        Lag(ix2).ContraMotorDist = MotorDiff(UniqueContraAxons(i),df);
        Lag(ix2).ContraSmallCount = smallCount;
        Lag(ix2).ContraLargeCount = largeCount;
        Lag(ix2).ContraLargeDiffSmall = largeCount-smallCount;
        ix2 = ix2+1;
    end
    clear A;
end

for i = 1:numel(DbxCells)
    tempLeadIDs = find(ismember(DBX(i).Inputs, [Lead.ContraAxonID]));
    tempLagIDs = find(ismember(DBX(i).Inputs, [Lag.ContraAxonID]));
    DBX(i).LeadContraPathLength = DBX(i).PathLength(tempLeadIDs)/max(Pvec_tree(DBX(i).Tree{1}));
    DBX(i).LeadContraPSD = DBX(i).PSDsize(tempLeadIDs);
    LeadContraPathLengths(i,:) = histcounts(DBX(i).LeadContraPathLength,0:0.1:1);
    
    DBX(i).LagContraPathLength = DBX(i).PathLength(tempLagIDs)/max(Pvec_tree(DBX(i).Tree{1}));
    DBX(i).LagContraPSD = DBX(i).PSDsize(tempLagIDs);
    LagContraPathLengths(i,:) = histcounts(DBX(i).LagContraPathLength,0:0.1:1);    
    
    clear tempLeadIDs;
    clear tempLagIDs;
end


subplot(4,4,13)
scatter([Lead.ContraDiff],[Lead.ContraMotorDist],20,'filled','MarkerFaceColor',leadColor);
hold on;
scatter([Lag.ContraDiff],[Lag.ContraMotorDist],20,'filled','MarkerFaceColor',lagColor);
line([-10,10],[0,0],'color','k','LineStyle','--');
line([0,0],[-100,100],'color','k','LineStyle','--');
axis square;
xlabel('Synapse diff');
ylabel('A-Ai diff');

subplot(4,4,14)
histogram(vertcat(DBX.LeadContraPathLength),10,'FaceColor',leadColor,'EdgeColor','none');
hold on;
histogram(vertcat(DBX.LagContraPathLength),10,'FaceColor',lagColor,'EdgeColor','none');
axis square;
box off;
xlabel('Norm. pathlength');
ylabel('count')

subplot(4,4,15)
shadedErrorBar(0.1:0.1:1,mean(LeadContraPathLengths),std(LeadContraPathLengths)/sqrt(9),'lineProps',{'Color',leadColor,'LineWidth',2});
hold on;
shadedErrorBar(0.1:0.1:1,mean(LagContraPathLengths),std(LagContraPathLengths)/sqrt(9),'lineProps',{'Color',lagColor,'LineWidth',2});
axis square;
xlabel('Norm pathlength');
ylabel('count per cell');

%% Everything Else module


UniqueEverythingElseAxons = unique(vertcat(DBX.EverythingElse));
    ix1 = 1;
    ix2 = 1;  
for i = 1:numel(UniqueEverythingElseAxons)
    [A,~] = SynapticPartners(UniqueEverythingElseAxons(i),2,df);
    %A = A(A<1e5);
    leadSum  = sum(ismember(A,leadNeurons));
    lagSum = sum(ismember(A,lagNeurons));
    smallCount = sum(ismember(A,smallABDCellIDs));
    largeCount = sum(ismember(A,largeABDCellIDs));

    if leadSum>lagSum
        Lead(ix1).EverythingElseAxonID = UniqueEverythingElseAxons(i);
        Lead(ix1).EverythingElseDiff = leadSum-lagSum;
        Lead(ix1).EverythingElseMotorDist = MotorDiff(UniqueEverythingElseAxons(i),df);
        Lead(ix1).EverythingElseSmallCount = smallCount;
        Lead(ix1).EverythingElseLargeCount = largeCount;
        Lead(ix1).EverythingElseLargeDiffSmall = largeCount-smallCount;
        ix1 = ix1+1;
    else
        Lag(ix2).EverythingElseAxonID = UniqueEverythingElseAxons(i);
        Lag(ix2).EverythingElseDiff = leadSum-lagSum;
        Lag(ix2).EverythingElseMotorDist = MotorDiff(UniqueEverythingElseAxons(i),df);
        Lag(ix2).EverythingElseSmallCount = smallCount;
        Lag(ix2).EverythingElseLargeCount = largeCount;
        Lag(ix2).EverythingElseLargeDiffSmall = largeCount-smallCount;
        ix2 = ix2+1;
    end
    clear A;
end

for i = 1:numel(DbxCells)
    tempLeadIDs = find(ismember(DBX(i).Inputs, [Lead.EverythingElseAxonID]));
    tempLagIDs = find(ismember(DBX(i).Inputs, [Lag.EverythingElseAxonID]));
    DBX(i).LeadEverythingElsePathLength = DBX(i).PathLength(tempLeadIDs)/max(Pvec_tree(DBX(i).Tree{1}));
    DBX(i).LeadEverythingElsePSD = DBX(i).PSDsize(tempLeadIDs);
    LeadEverythingElsePathLengths(i,:) = histcounts(DBX(i).LeadEverythingElsePathLength,0:0.1:1);
    
    DBX(i).LagEverythingElsePathLength = DBX(i).PathLength(tempLagIDs)/max(Pvec_tree(DBX(i).Tree{1}));
    DBX(i).LagEverythingElsePSD = DBX(i).PSDsize(tempLagIDs);
    LagEverythingElsePathLengths(i,:) = histcounts(DBX(i).LagEverythingElsePathLength,0:0.1:1);    
    
    clear tempLeadIDs;
    clear tempLagIDs;
end


% subplot(4,4,5)
% scatter([Lead.EverythingElseDiff],[Lead.EverythingElseMotorDist],20,'filled','MarkerFaceColor',leadColor);
% hold on;
% scatter([Lag.EverythingElseDiff],[Lag.EverythingElseMotorDist],20,'filled','MarkerFaceColor',lagColor);
% line([-10,10],[0,0],'color','k','LineStyle','--');
% line([0,0],[-100,100],'color','k','LineStyle','--');
% axis square;
% xlabel('Synapse diff');
% ylabel('A-Ai diff');
% 
% subplot(4,4,6)
% histogram(vertcat(DBX.LeadEverythingElsePathLength),10,'FaceColor',leadColor,'EdgeColor','none');
% hold on;
% histogram(vertcat(DBX.LagEverythingElsePathLength),10,'FaceColor',lagColor,'EdgeColor','none');
% axis square;
% box off;
% xlabel('Norm. pathlength');
% ylabel('count')
% 
% subplot(4,4,7)
% shadedErrorBar(0.1:0.1:1,mean(LeadEverythingElsePathLengths),std(LeadEverythingElsePathLengths)/sqrt(9),'lineProps',{'Color',leadColor,'LineWidth',2});
% hold on;
% shadedErrorBar(0.1:0.1:1,mean(LagEverythingElsePathLengths),std(LagEverythingElsePathLengths)/sqrt(9),'lineProps',{'Color',lagColor,'LineWidth',2});
% axis square;
% xlabel('Norm pathlength');
% ylabel('count per cell');


