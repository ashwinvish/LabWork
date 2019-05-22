clear;
addpath(genpath('/Users/ashwin/Documents/'));

colors = cbrewer('qual','Dark2',10);
startup

temp1 = colorcet('CBTL1','N',5); %  linear-tritanopic_krjcw_5-98_c46_n256
temp2 = colorcet('CBL2','N',5);  %  linear-protanopic-deuteranopic_kbw_5-98_c40_n256

leadColor = temp1(3,:); % dark red, CB safe
lagColor = temp2(3,:); % dark blue, CB safe

if ismac
    addpath(genpath('/Users/ashwin/Documents/'));
    df = readtable('/Users/ashwin/Documents/SynapseDetector/04152019.csv');
else
    addpath(genpath('/usr/people/ashwinv/seungmount/research/Ashwin/Scripts'));
    df = readtable('/usr/people/ashwinv/seungmount/research/Ashwin/SynapseDetector/11252018.csv');
end

confirmedALX = [76181 76201 76187 76184 76192 76197 ];
putativeALX = [76657 76623 76701 76661 76673 76690 76677 76624 77327 77580 ...
    77341 77344 77868 77872 77349 76634 77592 77578 77607 78629 77581 77591 ...
    77602 77821 77844 77845 77822 78406 78667 79022 79033 78404 78441 78421 ...
    79044 79046 79048 80221 78853 79017 79852 78451 79042 80596 80606 78911 ...
    79746 80271 79720 79976 77586 77369 78633 80750 77142 79060 78453 80885];

allALX = [confirmedALX,putativeALX];

for i = 1:numel(allALX)
    ALX(i) = InputsByClass(allALX(i),df);
end
%%
index =1;
for i = 1:numel(ALX)
    if ~isempty(ALX(i).Origin)
   alxOrigins(index,:) = ALX(i).Origin;
   index = index+1;
    end
end

%DBXpariwiseDist = pdist(dbxOrigins);
%[~,RCindex] = sort(DBXpariwiseDist);

for i = 1:size(alxOrigins,1)
    for j = 1:size(alxOrigins,1)
        ALXpariwiseDist(i,j) = sqrt(sum((alxOrigins(i,:) - alxOrigins(j,:)).^2));
        commonSaccadic = intersect(ALX(i).Saccadic,ALX(j).Saccadic);
        commonVestibular = intersect(ALX(i).Vestibular,ALX(j).Vestibular);
        commonIntgrators = intersect(ALX(i).Integrator,ALX(j).Integrator);
        commonContra = intersect(ALX(i).Contra,ALX(j).Contra);
        
        ALXcommonSaccadic(i,j) = length(commonSaccadic);
        ALXcommonVestibular(i,j) = length(commonVestibular);
        ALXcommonIntegrators(i,j) = length(commonIntgrators);
        ALXcommonContra(i,j) = length(commonContra);
        
        clear commonSaccadic;
        clear commonVestibular;
        clear commonIntgrators;
        clear commonContra;
    end
end

ALXcommonSaccadic(ALXcommonSaccadic ==0) = NaN;
ALXcommonVestibular(ALXcommonVestibular ==0) = NaN;
ALXcommonIntegrators(ALXcommonIntegrators ==0) = NaN;
ALXcommonContra(ALXcommonContra ==0) = NaN;
ALXpariwiseDist(ALXpariwiseDist ==0) = NaN;
[h,e,ind] = histcounts(ALXpariwiseDist(:),0:10:120);

figure;
subplot(4,4,1)
for i = 1:max(ind)
    errorbar(nanmean(ALXpariwiseDist(i==ind)),nanmean(ALXcommonSaccadic(i==ind)),...
        nanstd(ALXcommonSaccadic(i==ind))./sqrt(numel(ALXcommonSaccadic(i==ind))),'vertical',...
        'o','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',4,'Color','k','LineWidth',2);
    hold on
    
end

axis square
xlabel('Pairwise distance (\mum)');
ylabel('Common axons');
title('Saccadic')
box off;

subplot(4,4,2)
for i = 1:max(ind)
    errorbar(nanmean(ALXpariwiseDist(i==ind)),nanmean(ALXcommonVestibular(i==ind)),...
        nanstd(ALXcommonVestibular(i==ind))./sqrt(numel(ALXcommonVestibular(i==ind))),'vertical',...
        'o','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',4,'Color','k','LineWidth',2);
    hold on
    
end

axis square
xlabel('Pairwise distance (\mum)');
ylabel('Common axons');
title('Vestibular')
box off;

subplot(4,4,3)
for i = 1:max(ind)
    errorbar(nanmean(ALXpariwiseDist(i==ind)),nanmean(ALXcommonIntegrators(i==ind)),...
        nanstd(ALXcommonIntegrators(i==ind))./sqrt(numel(ALXcommonIntegrators(i==ind))),'vertical',...
        'o','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',4,'Color','k','LineWidth',2);
    hold on
    
end

axis square
xlabel('Pairwise distance (\mum)');
ylabel('Common axons');
title('Integrators')
box off;


subplot(4,4,4)
for i = 1:max(ind)
    errorbar(nanmean(ALXpariwiseDist(i==ind)),nanmean(ALXcommonContra(i==ind)),...
        nanstd(ALXcommonContra(i==ind))./sqrt(numel(ALXcommonContra(i==ind))),'vertical',...
        'o','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',4,'Color','k','LineWidth',2);
    hold on
    
end

axis square
xlabel('Pairwise distance (\mum)');
ylabel('Common axons');
title('Contra')
box off;

%%

uniqeALXSaccadicAxons = unique(vertcat(ALX.Saccadic));
ALXSaccadicMotorDist = MotorDiff(uniqeALXSaccadicAxons,df);

LeadLikeALXSaccadicAxons = uniqeALXSaccadicAxons(ALXSaccadicMotorDist>0);
LagLikeALXSaccadicAxons = uniqeALXSaccadicAxons(ALXSaccadicMotorDist<0);

save('LeadLikeALXSaccadicAxons.mat','LeadLikeALXSaccadicAxons');
save('LagLikeALXSaccadicAxons.mat','LagLikeALXSaccadicAxons');

LeadLikeALX = [];
LagLikeALX = [];
for i = 1:numel(allALX)
    [a,~] = SynapticPartners(allALX(i),1,df);
    a = a(a<1e5);
    leadSum = 0;
    lagSum = 0;
    for jj = 1:length(a)
        leadSum = leadSum + sum(a(jj) == LeadLikeALXSaccadicAxons);
        lagSum = lagSum + sum(a(jj) == LagLikeALXSaccadicAxons);
    end
    if (leadSum-lagSum) > 0 
        LeadLikeALX = [LeadLikeALX; allALX(i)];
    else
        LagLikeALX = [LagLikeALX; allALX(i)];
    end
end

save('LeadLikeALX.mat','LeadLikeALX');
save('LagLikeALX.mat','LagLikeALX');
%histogram(ALXSaccadicMotorDist,10);


% subplot(2,3,[1,4])
% transform_swc_AV(LeadLikeALXSaccadicAxons,colors(1,:),[],true,false);
% subplot(2,3,[2,5])
% transform_swc_AV(LagLikeALXSaccadicAxons,colors(2,:),[],true,false);
%% Saccadic Module
load('smallABDneurons.mat');
load('largeABDneurons.mat');

UniqueSaccadicAxons = unique(vertcat(ALX.Saccadic));
    ix1 = 1;
    ix2 = 1;
    
for i = 1:numel(UniqueSaccadicAxons)
    [A,~] = SynapticPartners(UniqueSaccadicAxons(i),2,df);
    A = A(A<1e5);
    leadSum  = sum(ismember(A,LeadLikeALX));
    lagSum = sum(ismember(A,LagLikeALX));
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

for i = 1:numel(ALX)
    tempLeadIDs = find(ismember(ALX(i).Inputs, [Lead.SaccadeAxonID]));
    tempLagIDs = find(ismember(ALX(i).Inputs, [Lag.SaccadeAxonID]));
    ALX(i).LeadSaccadicPathLength = ALX(i).PathLength(tempLeadIDs)/max(Pvec_tree(ALX(i).Tree{1}));
    ALX(i).LeadSaccadicPSD = ALX(i).PSDsize(tempLeadIDs);
    LeadSaccadicPathLengths(i,:) = histcounts(ALX(i).LeadSaccadicPathLength,0:0.1:1);
    
    ALX(i).LagSaccadicPathLength = ALX(i).PathLength(tempLagIDs)/max(Pvec_tree(ALX(i).Tree{1}));
    ALX(i).LagSaccadicPSD = ALX(i).PSDsize(tempLagIDs);
    LagSaccadicPathLengths(i,:) = histcounts(ALX(i).LagSaccadicPathLength,0:0.1:1);    
    
    clear tempLeadIDs;
    clear tempLagIDs;
end


figure;
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
histogram(vertcat(ALX.LeadSaccadicPathLength),10,'FaceColor',leadColor,'EdgeColor','none');
hold on;
histogram(vertcat(ALX.LagSaccadicPathLength),10,'FaceColor',lagColor,'EdgeColor','none');
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


UniqueVestibularAxons = unique(vertcat(ALX.Vestibular));
    ix1 = 1;
    ix2 = 1;  
for i = 1:numel(UniqueVestibularAxons)
    [A,~] = SynapticPartners(UniqueVestibularAxons(i),2,df);
    %A = A(A<1e5);
    leadSum  = sum(ismember(A,LeadLikeALX));
    lagSum = sum(ismember(A,LagLikeALX));
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

for i = 1:numel(ALX)
    tempLeadIDs = find(ismember(ALX(i).Inputs, [Lead.VestibularAxonID]));
    tempLagIDs = find(ismember(ALX(i).Inputs, [Lag.VestibularAxonID]));
    ALX(i).LeadVestibularPathLength = ALX(i).PathLength(tempLeadIDs)/max(Pvec_tree(ALX(i).Tree{1}));
    ALX(i).LeadVestibularPSD = ALX(i).PSDsize(tempLeadIDs);
    LeadVestibularPathLengths(i,:) = histcounts(ALX(i).LeadVestibularPathLength,0:0.1:1);
    
    ALX(i).LagVestibularPathLength = ALX(i).PathLength(tempLagIDs)/max(Pvec_tree(ALX(i).Tree{1}));
    ALX(i).LagVestibularPSD = ALX(i).PSDsize(tempLagIDs);
    LagVestibularPathLengths(i,:) = histcounts(ALX(i).LagVestibularPathLength,0:0.1:1);    
    
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
histogram(vertcat(ALX.LeadVestibularPathLength),10,'FaceColor',leadColor,'EdgeColor','none');
hold on;
histogram(vertcat(ALX.LagVestibularPathLength),10,'FaceColor',lagColor,'EdgeColor','none');
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


UniqueIntegratorAxons = unique(vertcat(ALX.Integrator));
    ix1 = 1;
    ix2 = 1;  
for i = 1:numel(UniqueIntegratorAxons)
    [A,~] = SynapticPartners(UniqueIntegratorAxons(i),2,df);
    %A = A(A<1e5);
    leadSum  = sum(ismember(A,LeadLikeALX));
    lagSum = sum(ismember(A,LagLikeALX));
    smallCount = sum(ismember(A,smallABDCellIDs));
    largeCount = sum(ismember(A,largeABDCellIDs));

    if leadSum>lagSum
        Lead(ix1).IntegratorAxonID = UniqueIntegratorAxons(i);
        Lead(ix1).IntegratorLeadSum = leadSum;
        Lead(ix1).IntegratorLagSum = lagSum;
        Lead(ix1).IntegratorDiff = leadSum-lagSum;
        Lead(ix1).IntegratorMotorDist = MotorDiff(UniqueIntegratorAxons(i),df);
        Lead(ix1).IntegratorSmallCount = smallCount;
        Lead(ix1).IntegratorLargeCount = largeCount;
        Lead(ix1).IntegratorLargeDiffSmall = largeCount-smallCount;
        ix1 = ix1+1;
    else
        Lag(ix2).IntegratorAxonID = UniqueIntegratorAxons(i);
        Lag(ix2).IntegratorLeadSum = leadSum;
        Lag(ix2).IntegratorLagSum = lagSum;
        Lag(ix2).IntegratorDiff = leadSum-lagSum;
        Lag(ix2).IntegratorMotorDist = MotorDiff(UniqueIntegratorAxons(i),df);
        Lag(ix2).IntegratorSmallCount = smallCount;
        Lag(ix2).IntegratorLargeCount = largeCount;
        Lag(ix2).IntegratorLargeDiffSmall = largeCount-smallCount;
        ix2 = ix2+1;
    end
    clear A;
end

for i = 1:numel(ALX)
    tempLeadIDs = find(ismember(ALX(i).Inputs, [Lead.IntegratorAxonID]));
    tempLagIDs = find(ismember(ALX(i).Inputs, [Lag.IntegratorAxonID]));
    ALX(i).LeadIntegratorPathLength = ALX(i).PathLength(tempLeadIDs)/max(Pvec_tree(ALX(i).Tree{1}));
    ALX(i).LeadIntegratorPSD = ALX(i).PSDsize(tempLeadIDs);
    LeadIntegratorPathLengths(i,:) = histcounts(ALX(i).LeadIntegratorPathLength,0:0.1:1);
    
    ALX(i).LagIntegratorPathLength = ALX(i).PathLength(tempLagIDs)/max(Pvec_tree(ALX(i).Tree{1}));
    ALX(i).LagIntegratorPSD = ALX(i).PSDsize(tempLagIDs);
    LagIntegratorPathLengths(i,:) = histcounts(ALX(i).LagIntegratorPathLength,0:0.1:1);    
    
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
histogram(vertcat(ALX.LeadIntegratorPathLength),10,'FaceColor',leadColor,'EdgeColor','none');
hold on;
histogram(vertcat(ALX.LagIntegratorPathLength),10,'FaceColor',lagColor,'EdgeColor','none');
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

subplot(4,4,12) 

plot(repmat([1;2],1,size([Lead.IntegratorLeadSum],2)),[[Lead.IntegratorLeadSum];[Lead.IntegratorLagSum]],'o','Color','w');
hold on
plot([1,2],[mean()]

violinPlot([[Lead.IntegratorLeadSum]',[Lead.IntegratorLagSum]'],'showMM',2,'color',leadColor);

plot(repmat([4;5],1,size([Lag.IntegratorLeadSum],2)),[[Lag.IntegratorLeadSum];[Lag.IntegratorLagSum]],'-o','Color',lagColor);
boxplot([Lag.IntegratorLeadSum],4*ones(1,size([Lag.IntegratorLeadSum],2)),'BoxStyle','filled','Colors',lagColor,'Jitter',2,'MedianStyle','target');
boxplot([Lag.IntegratorLeadSum],5*ones(1,size([Lag.IntegratorLagSum],2)),'BoxStyle','filled','Colors',lagColor,'Jitter',2,'MedianStyle','target');

set(gca,'XLim',[0,6]);
box off



%% Contra Module

UniqueContraAxons = unique(vertcat(ALX.Contra));
    ix1 = 1;
    ix2 = 1;  
for i = 1:numel(UniqueContraAxons)
    [A,~] = SynapticPartners(UniqueContraAxons(i),2,df);
    %A = A(A<1e5);
    leadSum  = sum(ismember(A,LeadLikeALX));
    lagSum = sum(ismember(A,LagLikeALX));
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

for i = 1:numel(ALX)
    tempLeadIDs = find(ismember(ALX(i).Inputs, [Lead.ContraAxonID]));
    tempLagIDs = find(ismember(ALX(i).Inputs, [Lag.ContraAxonID]));
    ALX(i).LeadContraPathLength = ALX(i).PathLength(tempLeadIDs)/max(Pvec_tree(ALX(i).Tree{1}));
    ALX(i).LeadContraPSD = ALX(i).PSDsize(tempLeadIDs);
    LeadContraPathLengths(i,:) = histcounts(ALX(i).LeadContraPathLength,0:0.1:1);
    
    ALX(i).LagContraPathLength = ALX(i).PathLength(tempLagIDs)/max(Pvec_tree(ALX(i).Tree{1}));
    ALX(i).LagContraPSD = ALX(i).PSDsize(tempLagIDs);
    LagContraPathLengths(i,:) = histcounts(ALX(i).LagContraPathLength,0:0.1:1);    
    
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
histogram(vertcat(ALX.LeadContraPathLength),10,'FaceColor',leadColor,'EdgeColor','none');
hold on;
histogram(vertcat(ALX.LagContraPathLength),10,'FaceColor',lagColor,'EdgeColor','none');
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


UniqueEverythingElseAxons = unique(vertcat(ALX.EverythingElse));
    ix1 = 1;
    ix2 = 1;  
for i = 1:numel(UniqueEverythingElseAxons)
    [A,~] = SynapticPartners(UniqueEverythingElseAxons(i),2,df);
    %A = A(A<1e5);
    leadSum  = sum(ismember(A,LeadLikeALX));
    lagSum = sum(ismember(A,LagLikeALX));
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

for i = 1:numel(ALX)
    tempLeadIDs = find(ismember(ALX(i).Inputs, [Lead.EverythingElseAxonID]));
    tempLagIDs = find(ismember(ALX(i).Inputs, [Lag.EverythingElseAxonID]));
    ALX(i).LeadEverythingElsePathLength = ALX(i).PathLength(tempLeadIDs)/max(Pvec_tree(ALX(i).Tree{1}));
    ALX(i).LeadEverythingElsePSD = ALX(i).PSDsize(tempLeadIDs);
    LeadEverythingElsePathLengths(i,:) = histcounts(ALX(i).LeadEverythingElsePathLength,0:0.1:1);
    
    ALX(i).LagEverythingElsePathLength = ALX(i).PathLength(tempLagIDs)/max(Pvec_tree(ALX(i).Tree{1}));
    ALX(i).LagEverythingElsePSD = ALX(i).PSDsize(tempLagIDs);
    LagEverythingElsePathLengths(i,:) = histcounts(ALX(i).LagEverythingElsePathLength,0:0.1:1);    
    
    clear tempLeadIDs;
    clear tempLagIDs;
end



%% partners of lead and Lag axons and how recurrently they.

[LeadLikeAxons,countsLead,Leadgraph,LeadShortestPaths] = PartnerConnectivity(LeadLikeALXSaccadicAxons,uniqeALXSaccadicAxons,'Saccadic',df);
[LagLikeAxons,countsLag,Laggraph,LagShortestPaths] = PartnerConnectivity(LagLikeALXSaccadicAxons,uniqeALXSaccadicAxons,'Saccadic',df);

countsLeadLag = vertcat(countsLead,countsLag);

subplot(2,3,3)
imagesc(countsLeadLag);
minMax = [min(countsLeadLag(:)), max(countsLeadLag(:))];
colorcet('L1','N',minMax(2)-minMax(1),'reverse',1);
line([0,size(countsLeadLag,2)],[size(LeadLikeAxons,2)+0.5 ,size(LeadLikeAxons,2)+0.5],'color','k','LineWidth',2);
line([size(LeadLikeAxons,2)+0.5 ,size(LeadLikeAxons,2)+0.5],[0,size(countsLeadLag,2)],'color','k','LineWidth',2);
colorbar;
box on;
daspect([1,1,1]);

subplot(2,3,6)
histogram(LeadShortestPaths,'FaceColor','none','EdgeColor',colors(1,:),'LineWidth',4,'DisplayStyle','stairs');
hold on;
histogram(LagShortestPaths,'FaceColor','none','EdgeColor',colors(2,:),'LineWidth',4,'DisplayStyle','stairs')
box off;
xlabel('shortest path nodes');
ylabel('count');

%%






