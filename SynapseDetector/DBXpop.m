% DBX population distributions
clear;
addpath(genpath('/Users/ashwin/Documents/'));

colors = cbrewer('qual','Dark2',10);

temp1 = colorcet('CBTL1','N',5); %  linear-tritanopic_krjcw_5-98_c46_n256
temp2 = colorcet('CBL2','N',5);  %  linear-protanopic-deuteranopic_kbw_5-98_c40_n256

leadColor = temp1(3,:); % dark red, CB safe
lagColor = temp2(3,:); % dark blue, CB safe

PartnerColors = colorcet('CBTD1','N',10);
lightRed = PartnerColors(10,:);
lightBlue = PartnerColors(1,:);


startup

if ismac
    addpath(genpath('/Users/ashwin/Documents/'));
    df = readtable('/Users/ashwin/Documents/SynapseDetector/11252018.csv');
else
    addpath(genpath('/usr/people/ashwinv/seungmount/research/Ashwin/Scripts'));
    df = readtable('/usr/people/ashwinv/seungmount/research/Ashwin/SynapseDetector/11252018.csv');
end


confirmedDBX = [76199 76200 76182 76183 76185 76186 76189 76191 76188 ];
putativeDBX = [77582 77588 77595 77597 77598 77599 77605 79040 76399 76828 ...
    76829 78435 76826 76874 76289 76542 76832 76836 76838 76877 76883 76884 ...
    78400 78429 76383 76867 77594 76830 76876 78440 78438 78447 79045 78434 ...
    78450 76878 77683 77450 78448 ];

AllDBX = [confirmedDBX, putativeDBX];

for i = 1:numel(AllDBX)
    DBX(i) = InputsByClass(AllDBX(i),df);
end

% 
for i = 1:numel(AllDBX)
    DBXSaccadic(i).PathLength = DBX(i).PathLength(DBX(i).isSaccadic)/max(Pvec_tree(DBX(i).Tree{1}));
    DBXVestibular(i).PathLength = DBX(i).PathLength(DBX(i).isVestibular)/max(Pvec_tree(DBX(i).Tree{1}));
    DBXContra(i).PathLength = DBX(i).PathLength(DBX(i).isContra)/max(Pvec_tree(DBX(i).Tree{1}));
    DBXIntegrator(i).PathLength = DBX(i).PathLength(DBX(i).isIntegrator)/max(Pvec_tree(DBX(i).Tree{1}));
    
    DBXSaccadic(i).PSDsize = DBX(i).PSDsize(DBX(i).isSaccadic);
    DBXVestibular(i).PSDsize = DBX(i).PSDsize(DBX(i).isVestibular);
    DBXContra(i).PSDsize = DBX(i).PSDsize(DBX(i).isContra);
    DBXIntegrator(i).PSDsize = DBX(i).PSDsize(DBX(i).isIntegrator)
end


subplot(1,3,1)
[~,plotOrder] = sortrows(vertcat(DBX.Origin),2);
for i = 1:numel(AllDBX)
    scatter(DBXSaccadic(plotOrder(i)).PathLength,DBX(plotOrder(i)).Origin(2)*ones(size(DBXSaccadic(plotOrder(i)).PathLength,1),1),...
        0.05*(DBXSaccadic(plotOrder(i)).PSDsize),colors(1,:));
    hold on;
    %if ~isnan(mean(DBXSaccadic(plotOrder(i)).PathLength))
    %    plot(mean(DBXSaccadic(plotOrder(i)).PathLength),DBX(plotOrder(i)).Origin(2),'-ro','markerSize',0.01*(median(DBXSaccadic(plotOrder(i)).PSDsize)));
    %end
    scatter(DBXVestibular(plotOrder(i)).PathLength,DBX(plotOrder(i)).Origin(2)*ones(size(DBXVestibular(plotOrder(i)).PathLength,1),1),...
        0.05*(DBXVestibular(plotOrder(i)).PSDsize),colors(2,:));
    scatter(DBXContra(plotOrder(i)).PathLength,DBX(plotOrder(i)).Origin(2)*ones(size(DBXContra(plotOrder(i)).PathLength,1),1),...
        0.05*(DBXContra(plotOrder(i)).PSDsize),colors(3,:));
     scatter(DBXIntegrator(plotOrder(i)).PathLength,DBX(plotOrder(i)).Origin(2)*ones(size(DBXIntegrator(plotOrder(i)).PathLength,1),1),...
        0.05*(DBXIntegrator(plotOrder(i)).PSDsize),colors(4,:));
end
%daspect([1,1,1]);
set(gca,'YDir','reverse');
ylabel('RC positon in \mu');
xlabel('Normalized path length');


% subplot(1,3,2);
% [~,plotOrder] = sortrows(vertcat(DBX.Origin),1);
% for i = 1:numel(AllDBX)
%     scatter(DBXSaccadic(plotOrder(i)).PathLength,DBX(plotOrder(i)).Origin(1)*ones(size(DBXSaccadic(plotOrder(i)).PathLength,1),1),...
%         0.05*(DBXSaccadic(plotOrder(i)).PSDsize),'k');
%     hold on;
%     %if ~isnan(mean(DBXSaccadic(plotOrder(i)).PathLength))
%     %    plot(mean(DBXSaccadic(plotOrder(i)).PathLength),DBX(plotOrder(i)).Origin(2),'-ro','markerSize',0.01*(median(DBXSaccadic(plotOrder(i)).PSDsize)));
%     %end
% end
% %daspect([1,1,1]);
% set(gca,'YDir','reverse');
% ylabel('ML positon in \mu');
% xlabel('Normalized path length');
% 
% subplot(1,3,3)
% [~,plotOrder] = sortrows(vertcat(DBX.Origin),3);
% for i = 1:numel(AllDBX)
%     scatter(DBXSaccadic(plotOrder(i)).PathLength,DBX(plotOrder(i)).Origin(3)*ones(size(DBXSaccadic(plotOrder(i)).PathLength,1),1),...
%         0.05*(DBXSaccadic(plotOrder(i)).PSDsize),'k');
%     hold on;
%     %if ~isnan(mean(DBXSaccadic(plotOrder(i)).PathLength))
%     %    plot(mean(DBXSaccadic(plotOrder(i)).PathLength),DBX(plotOrder(i)).Origin(2),'-ro','markerSize',0.01*(median(DBXSaccadic(plotOrder(i)).PSDsize)));
%     %end
% end
% %daspect([1,1,1]);
% set(gca,'YDir','reverse');
% ylabel('DV positon in \mu');
% xlabel('Normalized path length');

%% Saccadic Axons

UniqueDBXSaccadicAxons = unique(vertcat(DBX.Saccadic));
DBXSaccadicMotorDiff = MotorDiff(UniqueDBXSaccadicAxons,df);
figure;
%histogram(DBXSaccadicMotorDiff,10);

LeadLikeDBXSaccadicAxons = UniqueDBXSaccadicAxons(DBXSaccadicMotorDiff>0);
LagLikeDBXSaccadicAxons = UniqueDBXSaccadicAxons(DBXSaccadicMotorDiff<0);

save('LeadLikeDBXSaccadicAxons.mat','LeadLikeDBXSaccadicAxons');
save('LagLikeDBXSaccadicAxons.mat','LagLikeDBXSaccadicAxons');

leadAxonOrder = find(ismember(UniqueDBXSaccadicAxons,LeadLikeDBXSaccadicAxons));
lagAxonOrder = find(ismember(UniqueDBXSaccadicAxons,LagLikeDBXSaccadicAxons));


LeadLikeDBX = []
LagLikeDBX = [];
for i = 1:numel(AllDBX)
    [a,~] = SynapticPartners(AllDBX(i),1,df);
    a = a(a<1e5);
    leadSum = 0;
    lagSum = 0;
    for jj = 1:length(a)
        leadSum = leadSum + sum(a(jj) == LeadLikeDBXSaccadicAxons);
        lagSum = lagSum + sum(a(jj) == LagLikeDBXSaccadicAxons);
    end
    if (leadSum-lagSum) > 0 
        LeadLikeDBX = [LeadLikeDBX; AllDBX(i)];
    else
        LagLikeDBX = [LagLikeDBX; AllDBX(i)];
    end
end

save('LeadLikeDBX.mat','LeadLikeDBX');
save('LagLikeDBX.mat','LagLikeDBX');


subplot(2,3,[1,4])
transform_swc_AV(LeadLikeDBX,leadColor,[],true,false);
%subplot(2,3,[2,5])
transform_swc_AV(LagLikeDBX,lagColor,[],true,false);

subplot(2,3,[2,5])
LeadDBXOrder = find(ismember(AllDBX,LeadLikeDBX));
LagDBXOrder = find(ismember(AllDBX,LagLikeDBX));
LeadDBXOrigins = vertcat(DBX(LeadDBXOrder).Origin);
LagDBXOrigins = vertcat(DBX(LagDBXOrder).Origin);
histogram(LeadDBXOrigins(:,2),30,'EdgeColor',leadColor,'Orientation','horizontal','DisplayStyle','stairs','LineWidth',2);
hold on;
histogram(LagDBXOrigins(:,2),30,'EdgeColor',lagColor,'Orientation','horizontal','DisplayStyle','stairs','LineWidth',2);
set(gca,'YDir','reverse','YLim',[0,1280]);
daspect([1,30,1]);


%% Saccadic module

load('smallABDneurons.mat');
load('largeABDneurons.mat');

UniqueSaccadicAxons = unique(vertcat(DBX.Saccadic));
    ix1 = 1;
    ix2 = 1;
    
for i = 1:numel(UniqueSaccadicAxons)
    [A,~] = SynapticPartners(UniqueSaccadicAxons(i),2,df);
    A = A(A<1e5);
    leadSum  = sum(ismember(A,LeadLikeDBX));
    lagSum = sum(ismember(A,LagLikeDBX));
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

for i = 1:numel(DBX)
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
    leadSum  = sum(ismember(A,LeadLikeDBX));
    lagSum = sum(ismember(A,LagLikeDBX));
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

for i = 1:numel(DBX)
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
    leadSum  = sum(ismember(A,LeadLikeDBX));
    lagSum = sum(ismember(A,LagLikeDBX));
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

for i = 1:numel(DBX)
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
    leadSum  = sum(ismember(A,LeadLikeDBX));
    lagSum = sum(ismember(A,LagLikeDBX));
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

for i = 1:numel(DBX)
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
    leadSum  = sum(ismember(A,LeadLikeDBX));
    lagSum = sum(ismember(A,LagLikeDBX));
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

for i = 1:numel(DBX)
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



%%
[DBXSaccadicAxons,SaccadicCounts,SaccadicGraph,SaccadicShortestPaths] = PartnerConnectivity(UniqueDBXSaccadicAxons,UniqueDBXSaccadicAxons,'Saccadic',df);

SaccadicCounts = SaccadicCounts.*~eye(size(SaccadicCounts));

figure;

subplot(2,3,[1,4])
transform_swc_AV(LeadLikeDBXSaccadicAxons,lightRed,[],true,false);
subplot(2,3,[2,5])
transform_swc_AV(LagLikeDBXSaccadicAxons,lightBlue,[],true,false);

subplot(2,3,3)
imagesc(SaccadicCounts([leadAxonOrder;lagAxonOrder],[leadAxonOrder;lagAxonOrder]))
minMax = [min(SaccadicCounts(:)), max(SaccadicCounts(:))];
colorcet('L1','N',minMax(2)-minMax(1),'reverse',1);
line([0,size(SaccadicCounts,2)],[size(leadAxonOrder,1)+0.5 ,size(leadAxonOrder,1)+0.5],'color','k','LineWidth',2);
line([size(leadAxonOrder,1)+0.5 ,size(leadAxonOrder,1)+0.5],[0,size(SaccadicCounts,2)],'color','k','LineWidth',2);
colorbar;
box on;
daspect([1,1,1]);




% for i = 1:size(SaccadicShortestPaths,2)
%     ii = SaccadicShortestPaths(i).ij(1);
%     jj = SaccadicShortestPaths(i).ij(2);
%    distMat(ii,jj) = SaccadicShortestPaths(i).dist;
% end
% distMat = distMat.*~eye(size(distMat));



%% Integrator Axons

UniqueDBXIntegratorAxons = unique(vertcat(DBX.Integrator));
DBXIntegratorMotorDiff = MotorDiff(UniqueDBXIntegratorAxons,df);
histogram(DBXIntegratorMotorDiff,10);

LeadLikeIntegratorAxons = UniqueDBXIntegratorAxons(DBXIntegratorMotorDiff>0);
LagLikeIntegratorAxons = UniqueDBXIntegratorAxons(DBXIntegratorMotorDiff<0);

[DBXIntegratorAxons,IntegratorCounts,IntegratorGraph,IntegratorShortestPaths] = PartnerConnectivity(UniqueDBXIntegratorAxons,UniqueDBXIntegratorAxons,'Integrator',df);
IntegratorCounts = IntegratorCounts.*~eye(size(IntegratorCounts));


leadAxonOrder = find(ismember(UniqueDBXIntegratorAxons,LeadLikeIntegratorAxons));
lagAxonOrder = find(ismember(UniqueDBXIntegratorAxons,LagLikeIntegratorAxons));

figure;
subplot(2,3,[1,4])
transform_swc_AV(LeadLikeIntegratorAxons,lightRed,[],true,false);
subplot(2,3,[2,5])
transform_swc_AV(LagLikeIntegratorAxons,lightBlue,[],true,false);
subplot(2,3,3)
imagesc(IntegratorCounts([leadAxonOrder;lagAxonOrder],[leadAxonOrder;lagAxonOrder]));
minMax = [min(IntegratorCounts(:)), max(IntegratorCounts(:))];
colorcet('L1','N',minMax(2)-minMax(1),'reverse',1);
line([0,size(IntegratorCounts,2)],[size(leadAxonOrder,1)+0.5 ,size(leadAxonOrder,1)+0.5],'color','k','LineWidth',2);
line([size(leadAxonOrder,1)+0.5 ,size(leadAxonOrder,1)+0.5],[0,size(IntegratorCounts,2)],'color','k','LineWidth',2);
colorbar;
box on;
daspect([1,1,1]);

% IntegratorLoopSize = zeros(size(IntegratorCounts));
% 
% for i = 1:size(IntegratorCounts,1)
%     ii = IntegratorShortestPaths(i).ij(1);
%     jj = IntegratorShortestPaths(i).ij(2);
%    IntegratorLoopSize(ii,jj) = IntegratorShortestPaths(i).size;
% end
% %distMat = distMat.*~eye(size(distMat));
% 
% subplot(2,3,6)
% histogram(IntegratorShortestPaths,'FaceColor','none','EdgeColor',colors(1,:),'LineWidth',4,'DisplayStyle','stairs');
% hold on;
% histogram(IntegratorShortestPaths,'FaceColor','none','EdgeColor',colors(2,:),'LineWidth',4,'DisplayStyle','stairs')
% box off;
% xlabel('shortest path nodes');
% ylabel('count');

%% Vestibular Axons

UniqueDBXVestibularAxons = unique(vertcat(DBX.Vestibular));
DBXVestibularMotorDiff = MotorDiff(UniqueDBXVestibularAxons,df);
histogram(DBXVestibularMotorDiff,10);

LeadLikeVestibularAxons = UniqueDBXVestibularAxons(DBXVestibularMotorDiff>0);
LagLikeVestibularAxons = UniqueDBXVestibularAxons(DBXVestibularMotorDiff<0);

[DBXVestibularAxons,VestibularCounts,VestibularGraph,VestibularShortestPaths] = PartnerConnectivity(UniqueDBXVestibularAxons,UniqueDBXVestibularAxons,'Vestibular',df);
VestibularCounts = VestibularCounts.*~eye(size(VestibularCounts));


leadAxonOrder = find(ismember(UniqueDBXVestibularAxons,LeadLikeVestibularAxons));
lagAxonOrder = find(ismember(UniqueDBXVestibularAxons,LagLikeVestibularAxons));

figure;
subplot(2,3,[1,4])
transform_swc_AV(LeadLikeVestibularAxons,lightRed,[],true,false);
subplot(2,3,[2,5])
transform_swc_AV(LagLikeVestibularAxons,lightBlue,[],true,false);
subplot(2,3,3)
imagesc(VestibularCounts([leadAxonOrder;lagAxonOrder],[leadAxonOrder;lagAxonOrder]));
minMax = [min(VestibularCounts(:)), max(VestibularCounts(:))];
colorcet('L1','N',minMax(2)-minMax(1),'reverse',1);
line([0,size(VestibularCounts,2)],[size(leadAxonOrder,1)+0.5 ,size(leadAxonOrder,1)+0.5],'color','k','LineWidth',2);
line([size(leadAxonOrder,1)+0.5 ,size(leadAxonOrder,1)+0.5],[0,size(VestibularCounts,2)],'color','k','LineWidth',2);
colorbar;
box on;
daspect([1,1,1]);

% IntegratorLoopSize = zeros(size(IntegratorCounts));
% 
% for i = 1:size(IntegratorCounts,1)
%     ii = IntegratorShortestPaths(i).ij(1);
%     jj = IntegratorShortestPaths(i).ij(2);
%    IntegratorLoopSize(ii,jj) = IntegratorShortestPaths(i).size;
% end
% %distMat = distMat.*~eye(size(distMat));
% 
% subplot(2,3,6)
% histogram(IntegratorShortestPaths,'FaceColor','none','EdgeColor',colors(1,:),'LineWidth',4,'DisplayStyle','stairs');
% hold on;
% histogram(IntegratorShortestPaths,'FaceColor','none','EdgeColor',colors(2,:),'LineWidth',4,'DisplayStyle','stairs')
% box off;
% xlabel('shortest path nodes');
% ylabel('count');