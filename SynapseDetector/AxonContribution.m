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
%%  plot population data according to groups
% 2 --> Vglut_negative; blue
% 3 --> Vglut_positive; red

figure(1)
subplot(4,4,1)
shadedErrorBar(t,mean(Firing(:,CellOrder==2),2),std(Firing(:,CellOrder==2),[],2),'lineprops',{'b'});
hold on ;
shadedErrorBar(t,mean(Firing(:,CellOrder==3),2),std(Firing(:,CellOrder==3),[],2),'lineprops',{'r'});
xlabel('Time (sec)');
ylabel('df/f');
axis square;
title('Clusters from type');

%% Coorelation analyis of the cells.
CorrVals = corr(Firing);
CorrVals = CorrVals - eye(size(CorrVals));
% cutoff = 0.5;
% Z = linkage(CorrVals,'complete','correlation')

subplot(4,4,2);
molOrder = [find(CellOrder ==1)'; find(CellOrder ==2)'; find(CellOrder ==3)'];
imagesc(CorrVals(molOrder,molOrder));
colormap(colorcet('L1'));
% set(gca,'XTick',1:size(order,1),'XTickLabel',{CellOrder(order)},'YTick',1:size(order,1),'YTickLabel',{CellOrder(order)});
colorbar;
axis square;



[idx,c] = kmeans(CorrVals',2,'Distance','correlation','MaxIter',100);
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
colormap(colorcet('L1'));
% set(gca,'XTick',1:size(order,1),'XTickLabel',{CellOrder(order)},'YTick',1:size(order,1),'YTickLabel',{CellOrder(order)});
colorbar;
axis square;
% title('Sorted by coorelations');

subplot(4,4,6);
shadedErrorBar(t,mean(Firing(:,leadCorrIDs),2),std(Firing(:,leadCorrIDs),[],2),'lineprops',{leadColor});
hold on;
shadedErrorBar(t,mean(Firing(:,lagCorrIDs),2),std(Firing(:,lagCorrIDs),[],2),'lineprops',{lagColor});
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
shadedErrorBar(t,mean(Firing(:,randomPop1),2),std(Firing(:,randomPop1),[],2),'lineprops',{leadColor});
hold on;
shadedErrorBar(t,mean(Firing(:,randomPop2),2),std(Firing(:,randomPop2),[],2),'lineprops',{lagColor});
axis square;
xlabel('Time (sec)');
ylabel('df/f');
title('Random clusters');

%% Extract Lead and Lag neuron partners from coorelation analysis.

% lead Neurons
for i = 1:numel(leadNeurons)
    Lead(i) = InputsByClass(leadNeurons(i),df);
end

% lag neurons
for i = 1:numel(lagNeurons)
    Lag(i)  = InputsByClass(lagNeurons(i),df);
end

% plot location on cells based in rhombomered organization.

for i = 1:numel(leadNeurons)
    Lead(i).r3(:,1) =  Lead(i).InputsRhombomeres(Lead(i).InputsRhombomeres(:,5) == 1,1);
    [~,Lead(i).r3(:,2)] = ismember(Lead(i).r3(:,1),Lead(i).Inputs);
    
    Lead(i).r4(:,1) =  Lead(i).InputsRhombomeres(Lead(i).InputsRhombomeres(:,6) == 1,1);
    [~,Lead(i).r4(:,2)] = ismember(Lead(i).r4(:,1),Lead(i).Inputs);
    
    Lead(i).r5(:,1) =  Lead(i).InputsRhombomeres(Lead(i).InputsRhombomeres(:,7) == 1,1);
    [~,Lead(i).r5(:,2)] = ismember(Lead(i).r5(:,1),Lead(i).Inputs);
    
    Lead(i).r6(:,1) =  Lead(i).InputsRhombomeres(Lead(i).InputsRhombomeres(:,8) == 1,1);
    [~,Lead(i).r6(:,2)] = ismember(Lead(i).r6(:,1),Lead(i).Inputs);
    
    Lead(i).r7(:,1) =  Lead(i).InputsRhombomeres(Lead(i).InputsRhombomeres(:,9) == 1,1);
    [~,Lead(i).r7(:,2)] = ismember(Lead(i).r7(:,1),Lead(i).Inputs);
    
end

for i = 1:numel(lagNeurons)
    Lag(i).r3(:,1) =  Lag(i).InputsRhombomeres(Lag(i).InputsRhombomeres(:,5) == 1,1);
    [~,Lag(i).r3(:,2)] = ismember(Lag(i).r3(:,1),Lag(i).Inputs);
    
    Lag(i).r4(:,1) =  Lag(i).InputsRhombomeres(Lag(i).InputsRhombomeres(:,6) == 1,1);
    [~,Lag(i).r4(:,2)] = ismember(Lag(i).r4(:,1),Lag(i).Inputs);
    
    Lag(i).r5(:,1) =  Lag(i).InputsRhombomeres(Lag(i).InputsRhombomeres(:,7) == 1,1);
    [~,Lag(i).r5(:,2)] = ismember(Lag(i).r5(:,1),Lag(i).Inputs);
    
    Lag(i).r6(:,1) =  Lag(i).InputsRhombomeres(Lag(i).InputsRhombomeres(:,8) == 1,1);
    [~,Lag(i).r6(:,2)] = ismember(Lag(i).r6(:,1),Lag(i).Inputs);
    
    Lag(i).r7(:,1) =  Lag(i).InputsRhombomeres(Lag(i).InputsRhombomeres(:,9) == 1,1);
    [~,Lag(i).r7(:,2)] = ismember(Lag(i).r7(:,1),Lag(i).Inputs);
    
end

figure(2)

subplot(4,4,1)

histogram(vertcat(Lead.PSDsize)/sum(vertcat(Lead.PSDsize)),'FaceColor',leadColor,'EdgeColor','none');
hold on;
histogram(vertcat(Lag.PSDsize)/sum(vertcat(Lag.PSDsize)),'FaceColor',lagColor,'EdgeColor','none');
box off;
axis square;
title('PSDsize/sum(PSDsize)');


leadAxons = vertcat(Lead.Inputs);
numberOfLeadAxons = numel(leadAxons);
leadAxons = leadAxons(leadAxons<1e5);

lagAxons = vertcat(Lag.Inputs);
numberOfLagAxons = numel(lagAxons);
lagAxons = lagAxons(lagAxons<1e5);

subplot(4,4,2)
histogram(histcounts(leadAxons,unique(leadAxons)),'FaceColor',leadColor,'EdgeColor','none');
hold on;
histogram(histcounts(lagAxons,unique(lagAxons)),'FaceColor',lagColor,'EdgeColor','none');
box off;
axis square;
title('Number of synapses');
%% Extract lead/lag axons identity


% lead/lag saccadic neurons onto Motor neurons

for i = 1:numel(leadNeurons)
    Lead(i).SaccadicMotorDist = isMotor(Lead(i).Inputs(Lead(i).isSaccadic,:),df);
    Lead(i).VestibularMotorDist = isMotor(Lead(i).Inputs(Lead(i).isVestibular,:),df);
    Lead(i).ContraMotorDist  = isMotor(Lead(i).Inputs(Lead(i).isContra,:),df);
    Lead(i).IntegratorMotorDist = isMotor(Lead(i).Inputs(Lead(i).isIntegrator,:),df);
    Lead(i).EverythingElseMotorDist = isMotor(Lead(i).Inputs(Lead(i).isEverythingElse,:),df);
end


for i = 1:numel(lagNeurons)
    Lag(i).SaccadicMotorDist = isMotor(Lag(i).Inputs(Lag(i).isSaccadic,:),df);
    Lag(i).VestibularMotorDist = isMotor(Lag(i).Inputs(Lag(i).isVestibular,:),df);
    Lag(i).ContraMotorDist = isMotor(Lag(i).Inputs(Lag(i).isContra,:),df);
    Lag(i).IntegratorMotorDist = isMotor(Lag(i).Inputs(Lag(i).isIntegrator,:),df);
    Lag(i).EverythingElseMotorDist = isMotor(Lag(i).Inputs(Lag(i).isEverythingElse,:),df);
end

subplot(4,4,3)
plot([1,2,3,4,5],[size(vertcat(Lead.Vestibular),1),size(vertcat(Lead.Integrator),1),size(vertcat(Lead.Saccadic),1), ...
    size(vertcat(Lead.Contra),1), size(vertcat(Lead.EverythingElse),1)],'-ro');
hold on;
plot([1,2,3,4,5],[size(vertcat(Lag.Vestibular),1), size(vertcat(Lag.Integrator),1),size(vertcat(Lag.Saccadic),1), ...
    size(vertcat(Lag.Contra),1),size(vertcat(Lag.EverythingElse),1)],'-bo');
box off;
axis square;
xticks([1,2,3,4,5]);
xticklabels({'Vestibular','Integrator','Saccadic','Contra','Remaining'});
xtickangle(45);
ylabel('Synapses');


%% Contribution to each cell

% Saccadic contribution

leadSaccadicAxons = vertcat(Lead.Saccadic);
lagSaccadicAxons = vertcat(Lag.Saccadic);
uniqueSaccadicAxons = unique([leadSaccadicAxons;lagSaccadicAxons]);
saccadicMotorDist = [vertcat(Lead.SaccadicMotorDist);vertcat(Lag.SaccadicMotorDist)];

figure;

for i = 1:numel(uniqueSaccadicAxons)
    if  sum(leadSaccadicAxons == uniqueSaccadicAxons(i)) - sum(lagSaccadicAxons == uniqueSaccadicAxons(i)) > 0
        subplot(5,4,1)
        l = find(saccadicMotorDist(:,1) == uniqueSaccadicAxons(i));
        leadDiff(i).Saccadic =sum(leadSaccadicAxons == uniqueSaccadicAxons(i)) - sum(lagSaccadicAxons == uniqueSaccadicAxons(i));
        leadMotorSynapseDiff(i).Saccadic = (saccadicMotorDist(l(1),2)+saccadicMotorDist(l(1),3)) - (saccadicMotorDist(l(1),4)+saccadicMotorDist(l(1),4));
        temp  = isPostSynapseMotor(uniqueSaccadicAxons(i),df);
        leadMotorDiff(i).Saccadic = (temp(1,2)+temp(1,3)) - (temp(1,4)+temp(1,5));
        leadDiffAxons(i).Saccadic = uniqueSaccadicAxons(i);
        plot([1,2],[leadDiff(i).Saccadic,leadMotorSynapseDiff(i).Saccadic],'-','Color',[1,0,0,0.1]);
        hold on;
    else
        subplot(5,4,1)
        lagDiff(i).Saccadic = sum(leadSaccadicAxons == uniqueSaccadicAxons(i)) - sum(lagSaccadicAxons == uniqueSaccadicAxons(i));
        l = find(saccadicMotorDist(:,1) == uniqueSaccadicAxons(i));
        lagMotorSynapseDiff(i).Saccadic = (saccadicMotorDist(l(1),2)+saccadicMotorDist(l(1),3)) - (saccadicMotorDist(l(1),4)+saccadicMotorDist(l(1),4));
        temp  = isPostSynapseMotor(uniqueSaccadicAxons(i),df);
        lagMotorDiff(i).Saccadic = (temp(1,2)+temp(1,3)) - (temp(1,4)+temp(1,5));
        lagDiffAxons(i).Saccadic =  uniqueSaccadicAxons(i);
        plot([1,2],[lagDiff(i).Saccadic,lagMotorSynapseDiff(i).Saccadic],'-','Color',[0,0,1,0.1]);
        hold on;
    end
end


plot([1,2],[mean([leadDiff.Saccadic]), mean([leadMotorSynapseDiff.Saccadic])],'-o','Color',[1,0,0]);
plot([1,2],[mean([lagDiff.Saccadic]), mean([lagMotorSynapseDiff.Saccadic])],'-o','Color',[0,0,1]);
set(gca,'XLim',[0,3]);
title('Saccadic axons');
box off;
axis square;

subplot(5,4,2)

scatter([leadDiff.Saccadic],[leadMotorSynapseDiff.Saccadic],20,'o','MarkerFaceColor',colors(1,:),'MarkerEdgeColor','none');
hold on;
scatter([lagDiff.Saccadic],[lagMotorSynapseDiff.Saccadic],20,'o','MarkerFaceColor',colors(1,:),'MarkerEdgeColor','none');
line([0,0],[-100,100],'linestyle',':','color','k');
line([-10,10],[0,0],'linestyle',':','color','k');
%ylabel('#(ABD - ABDi) synapses');

box off;
axis square;

subplot(5,4,3)

scatter([leadDiff.Saccadic],[leadMotorDiff.Saccadic],20,'o','MarkerFaceColor',colors(1,:),'MarkerEdgeColor','none');
hold on;
scatter([lagDiff.Saccadic],[lagMotorDiff.Saccadic],20,'o','MarkerFaceColor',colors(1,:),'MarkerEdgeColor','none');
line([0,0],[-20,20],'linestyle',':','color','k');
line([-10,10],[0,0],'linestyle',':','color','k');
%ylabel('#(ABD - ABDi) neurons');

box off;
axis square;

subplot(5,4,4)

scatter([leadMotorSynapseDiff.Saccadic],[leadMotorDiff.Saccadic],20,'o','MarkerFaceColor',leadColor,'MarkerEdgeColor','none');
hold on
scatter([lagMotorSynapseDiff.Saccadic],[lagMotorDiff.Saccadic],20,'o','MarkerFaceColor',lagColor,'MarkerEdgeColor','none');
line([0,0],[-20,20],'linestyle',':','color','k');
line([-200,200],[0,0],'linestyle',':','color','k');
%ylabel('#(ABD - ABDi) neurons');
%xlabel('#(ABD - ABDi) synapses')
text(150,-15,sprintf('n=%d',size(uniqueSaccadicAxons,1)));
box off;
axis square;


% Vestibular contribution


leadVestibularAxons = vertcat(Lead.Vestibular);
lagVestibularAxons = vertcat(Lag.Vestibular);
uniqueVestibularAxons = unique([leadVestibularAxons;lagVestibularAxons]);
vestibularMotorDist = [vertcat(Lead.VestibularMotorDist);vertcat(Lag.VestibularMotorDist)];


for i = 1:numel(uniqueVestibularAxons)
    if  sum(leadVestibularAxons == uniqueVestibularAxons(i)) - sum(lagVestibularAxons == uniqueVestibularAxons(i)) > 0
        subplot(5,4,5)
        l = find(vestibularMotorDist(:,1) == uniqueVestibularAxons(i));
        leadDiff(i).Vestibular =sum(leadVestibularAxons == uniqueVestibularAxons(i)) - sum(lagVestibularAxons == uniqueVestibularAxons(i));
        leadMotorSynapseDiff(i).Vestibular = (vestibularMotorDist(l(1),2)+vestibularMotorDist(l(1),3)) - (vestibularMotorDist(l(1),4)+vestibularMotorDist(l(1),4));
        temp  = isPostSynapseMotor(uniqueVestibularAxons(i),df);
        leadMotorDiff(i).Vestibular = (temp(1,2)+temp(1,3)) - (temp(1,4)+temp(1,5));
        leadDiffAxons(i).Vestibular = uniqueVestibularAxons(i);
        plot([1,2],[leadDiff(i).Vestibular,leadMotorSynapseDiff(i).Vestibular],'-','Color',[1,0,0,0.1]);
        hold on;
    else
        subplot(5,4,5)
        lagDiff(i).Vestibular = sum(leadVestibularAxons == uniqueVestibularAxons(i)) - sum(lagVestibularAxons == uniqueVestibularAxons(i));
        l = find(vestibularMotorDist(:,1) == uniqueVestibularAxons(i));
        lagMotorSynapseDiff(i).Vestibular = (vestibularMotorDist(l(1),2)+vestibularMotorDist(l(1),3)) - (vestibularMotorDist(l(1),4)+vestibularMotorDist(l(1),4));
         temp  = isPostSynapseMotor(uniqueVestibularAxons(i),df);
        lagMotorDiff(i).Vestibular = (temp(1,2)+temp(1,3)) - (temp(1,4)+temp(1,5));
        lagDiffAxons(i).Vestibular = uniqueVestibularAxons(i);
        plot([1,2],[lagDiff(i).Vestibular,lagMotorSynapseDiff(i).Vestibular],'-','Color',[0,0,1,0.1]);
        hold on;
    end
end

plot([1,2],[mean([leadDiff.Vestibular]), mean([leadMotorSynapseDiff.Vestibular])],'-o','Color',[1,0,0]);
plot([1,2],[mean([lagDiff.Vestibular]), mean([lagMotorSynapseDiff.Vestibular])],'-o','Color',[0,0,1]);
set(gca,'XLim',[0,3]);
title('Vestibular axons');
box off;
axis square;

subplot(5,4,6)
scatter([leadDiff.Vestibular],[leadMotorSynapseDiff.Vestibular],20,'o','MarkerFaceColor',colors(2,:),'MarkerEdgeColor','none');
hold on;
scatter([lagDiff.Vestibular],[lagMotorSynapseDiff.Vestibular],20,'o','MarkerFaceColor',colors(2,:),'MarkerEdgeColor','none');
line([0,0],[-100,100],'linestyle',':','color','k');
line([-10,10],[0,0],'linestyle',':','color','k');
%ylabel('#(ABD - ABDi) synapses');
box off;
axis square;

subplot(5,4,7)

scatter([leadDiff.Vestibular],[leadMotorDiff.Vestibular],20,'o','MarkerFaceColor',colors(2,:),'MarkerEdgeColor','none');
hold on;
scatter([lagDiff.Vestibular],[lagMotorDiff.Vestibular],20,'o','MarkerFaceColor',colors(2,:),'MarkerEdgeColor','none');
line([0,0],[-20,20],'linestyle',':','color','k');
line([-10,10],[0,0],'linestyle',':','color','k');
%ylabel('#(ABD - ABDi) neurons');

box off;
axis square;

subplot(5,4,8)

scatter([leadMotorSynapseDiff.Vestibular],[leadMotorDiff.Vestibular],20,'o','MarkerFaceColor',leadColor,'MarkerEdgeColor','none');
hold on
scatter([lagMotorSynapseDiff.Vestibular],[lagMotorDiff.Vestibular],20,'o','MarkerFaceColor',lagColor,'MarkerEdgeColor','none');
line([0,0],[-20,20],'linestyle',':','color','k');
line([-200,200],[0,0],'linestyle',':','color','k');
%ylabel('#(ABD - ABDi) neurons');
%xlabel('#(ABD - ABDi) synapses')
text(150,-15,sprintf('n=%d',size(uniqueVestibularAxons,1)));
box off;
axis square;

% Integrator contribution

leadIntegratorAxons = vertcat(Lead.Integrator);
lagIntegratorAxons = vertcat(Lag.Integrator);
uniqueIntegratorAxons = unique([leadIntegratorAxons;lagIntegratorAxons]);
integratorMotorDist = [vertcat(Lead.IntegratorMotorDist);vertcat(Lag.IntegratorMotorDist)];


for i = 1:numel(uniqueIntegratorAxons)
    if  sum(leadIntegratorAxons == uniqueIntegratorAxons(i)) - sum(lagIntegratorAxons == uniqueIntegratorAxons(i)) > 0
        subplot(5,4,9)
        l = find(integratorMotorDist(:,1) == uniqueIntegratorAxons(i));
        leadDiff(i).Integrator =sum(leadIntegratorAxons == uniqueIntegratorAxons(i)) - sum(lagIntegratorAxons == uniqueIntegratorAxons(i));
        leadMotorSynapseDiff(i).Integrator = (integratorMotorDist(l(1),2)+integratorMotorDist(l(1),3)) - (integratorMotorDist(l(1),4)+integratorMotorDist(l(1),4));
        temp  = isPostSynapseMotor(uniqueIntegratorAxons(i),df);
        leadMotorDiff(i).Integrator = (temp(1,2)+temp(1,3)) - (temp(1,4)+temp(1,5));
        leadDiffAxons(i).Integrator = uniqueIntegratorAxons(i);
        plot([1,2],[leadDiff(i).Integrator,leadMotorSynapseDiff(i).Integrator],'-','Color',[1,0,0,0.1]);
        hold on;
    else
        subplot(5,4,9)
        lagDiff(i).Integrator = sum(leadIntegratorAxons == uniqueIntegratorAxons(i)) - sum(lagIntegratorAxons == uniqueIntegratorAxons(i));
        l = find(integratorMotorDist(:,1) == uniqueIntegratorAxons(i));
        lagMotorSynapseDiff(i).Integrator = (integratorMotorDist(l(1),2)+integratorMotorDist(l(1),3)) - (integratorMotorDist(l(1),4)+integratorMotorDist(l(1),4));
        temp  = isPostSynapseMotor(uniqueIntegratorAxons(i),df);
        lagMotorDiff(i).Integrator = (temp(1,2)+temp(1,3)) - (temp(1,4)+temp(1,5));
        lagDiffAxons(i).Integrator = uniqueIntegratorAxons(i);
        plot([1,2],[lagDiff(i).Integrator,lagMotorSynapseDiff(i).Integrator],'-','Color',[0,0,1,0.1]);
        hold on;
    end
end

plot([1,2],[mean([leadDiff.Integrator]), mean([leadMotorSynapseDiff.Integrator])],'-o','Color',[1,0,0]);
plot([1,2],[mean([lagDiff.Integrator]), mean([lagMotorSynapseDiff.Integrator])],'-o','Color',[0,0,1]);
set(gca,'XLim',[0,3]);
title('Putative integrator axons');
box off;
axis square;


subplot(5,4,10)
scatter([leadDiff.Integrator],[leadMotorSynapseDiff.Integrator],20,'o','MarkerFaceColor',colors(3,:),'MarkerEdgeColor','none');
hold on;
scatter([lagDiff.Integrator],[lagMotorSynapseDiff.Integrator],20,'o','MarkerFaceColor',colors(3,:),'MarkerEdgeColor','none');
line([0,0],[-100,100],'linestyle',':','color','k');
line([-10,10],[0,0],'linestyle',':','color','k');
%ylabel('#(ABD - ABDi) synapses');
box off;
axis square;

subplot(5,4,11)

scatter([leadDiff.Integrator],[leadMotorDiff.Integrator],20,'o','MarkerFaceColor',colors(3,:),'MarkerEdgeColor','none');
hold on;
scatter([lagDiff.Integrator],[lagMotorDiff.Integrator],20,'o','MarkerFaceColor',colors(3,:),'MarkerEdgeColor','none');
line([0,0],[-20,20],'linestyle',':','color','k');
line([-10,10],[0,0],'linestyle',':','color','k');
%ylabel('#(ABD - ABDi) neurons');

box off;
axis square;

subplot(5,4,12)
scatter([leadMotorSynapseDiff.Integrator],[leadMotorDiff.Integrator],20,'o','MarkerFaceColor',leadColor,'MarkerEdgeColor','none');
hold on
scatter([lagMotorSynapseDiff.Integrator],[lagMotorDiff.Integrator],20,'o','MarkerFaceColor',lagColor,'MarkerEdgeColor','none');
line([0,0],[-20,20],'linestyle',':','color','k');
line([-200,200],[0,0],'linestyle',':','color','k');
%ylabel('#(ABD - ABDi) neurons');
%xlabel('#(ABD - ABDi) synapses');
text(150,-15,sprintf('n=%d',size(uniqueIntegratorAxons,1)));

box off;
axis square;


% Contra contribution

leadContraAxons = vertcat(Lead.Contra);
lagContraAxons = vertcat(Lag.Contra);
uniqueContraAxons = unique([leadContraAxons;lagContraAxons]);
contraMotorDist = [vertcat(Lead.ContraMotorDist);vertcat(Lag.ContraMotorDist)];

for i = 1:numel(uniqueContraAxons)
    if  sum(leadContraAxons == uniqueContraAxons(i)) - sum(lagContraAxons == uniqueContraAxons(i)) > 0
        subplot(5,4,13)
        l = find(contraMotorDist(:,1) == uniqueContraAxons(i));
        leadDiff(i).Contra =sum(leadContraAxons == uniqueContraAxons(i)) - sum(lagContraAxons == uniqueContraAxons(i));
        leadMotorSynapseDiff(i).Contra = (contraMotorDist(l(1),2)+contraMotorDist(l(1),3)) - (contraMotorDist(l(1),4)+contraMotorDist(l(1),4));
        temp  = isPostSynapseMotor(uniqueContraAxons(i),df);
        leadMotorDiff(i).Contra = (temp(1,2)+temp(1,3)) - (temp(1,4)+temp(1,5));
        leadDiffAxons(i).Contra = uniqueContraAxons(i);
        plot([1,2],[leadDiff(i).Contra,leadMotorSynapseDiff(i).Contra],'-','Color',[1,0,0,0.1]);
        hold on;
    else
        subplot(5,4,13)
        lagDiff(i).Contra = sum(leadContraAxons == uniqueContraAxons(i)) - sum(lagContraAxons == uniqueContraAxons(i));
        l = find(contraMotorDist(:,1) == uniqueContraAxons(i));
        lagMotorSynapseDiff(i).Contra = (contraMotorDist(l(1),2)+contraMotorDist(l(1),3)) - (contraMotorDist(l(1),4)+contraMotorDist(l(1),4));
        temp  = isPostSynapseMotor(uniqueContraAxons(i),df);
        lagMotorDiff(i).Contra = (temp(1,2)+temp(1,3)) - (temp(1,4)+temp(1,5));
        lagDiffAxons(i).Contra = uniqueContraAxons(i);
        plot([1,2],[lagDiff(i).Contra,lagMotorSynapseDiff(i).Contra],'-','Color',[0,0,1,0.1]);
        hold on;
    end
end


plot([1,2],[mean([leadDiff.Contra]), mean([leadMotorSynapseDiff.Contra])],'-o','Color',[1,0,0]);
plot([1,2],[mean([lagDiff.Contra]), mean([lagMotorSynapseDiff.Contra])],'-o','Color',[0,0,1]);
set(gca,'XLim',[0,3]);
title('Contra axons');
box off;
axis square;


subplot(5,4,14)
scatter([leadDiff.Contra],[leadMotorSynapseDiff.Contra],20,'o','MarkerFaceColor',colors(4,:),'MarkerEdgeColor','none');
hold on;
scatter([lagDiff.Contra],[lagMotorSynapseDiff.Contra],20,'o','MarkerFaceColor',colors(4,:),'MarkerEdgeColor','none');
line([0,0],[-100,100],'linestyle',':','color','k');
line([-10,10],[0,0],'linestyle',':','color','k');
%ylabel('#(ABD - ABDi) synapses');

box off;
axis square;

subplot(5,4,15)

scatter([leadDiff.Contra],[leadMotorDiff.Contra],20,'o','MarkerFaceColor',colors(4,:),'MarkerEdgeColor','none');
hold on;
scatter([lagDiff.Contra],[lagMotorDiff.Contra],20,'o','MarkerFaceColor',colors(4,:),'MarkerEdgeColor','none');
line([0,0],[-20,20],'linestyle',':','color','k');
line([-10,10],[0,0],'linestyle',':','color','k');
%ylabel('#(ABD - ABDi) neurons');

box off;
axis square;

subplot(5,4,16)
scatter([leadMotorSynapseDiff.Contra],[leadMotorDiff.Contra],20,'o','MarkerFaceColor',leadColor,'MarkerEdgeColor','none');
hold on
scatter([lagMotorSynapseDiff.Contra],[lagMotorDiff.Contra],20,'o','MarkerFaceColor',lagColor,'MarkerEdgeColor','none');
line([0,0],[-20,20],'linestyle',':','color','k');
line([-200,200],[0,0],'linestyle',':','color','k');
%ylabel('#(ABD - ABDi) neurons');
%xlabel('#(ABD - ABDi) synapses');
text(150,-15,sprintf('n=%d',size(uniqueContraAxons,1)));

box off;
axis square;

% EverythginElse contribution

leadEverythingElseAxons = vertcat(Lead.EverythingElse);
lagEverythingElseAxons = vertcat(Lag.EverythingElse);
uniqueEverythingElseAxons = unique([leadEverythingElseAxons;lagEverythingElseAxons]);
everythginElseMotorDist = [vertcat(Lead.EverythingElseMotorDist);vertcat(Lag.EverythingElseMotorDist)];

for i = 1:numel(uniqueEverythingElseAxons)
    if  sum(leadEverythingElseAxons == uniqueEverythingElseAxons(i)) - sum(lagEverythingElseAxons == uniqueEverythingElseAxons(i)) > 0
        subplot(5,4,17)
        l = find(everythginElseMotorDist(:,1) ==uniqueEverythingElseAxons(i));
        leadDiff(i).EverythingElse =sum(leadEverythingElseAxons == uniqueEverythingElseAxons(i)) - sum(lagEverythingElseAxons == uniqueEverythingElseAxons(i));
        leadMotorSynapseDiff(i).EverythingElse = (everythginElseMotorDist(l(1),2)+everythginElseMotorDist(l(1),3)) - (everythginElseMotorDist(l(1),4)+everythginElseMotorDist(l(1),4));
        temp  = isPostSynapseMotor(uniqueEverythingElseAxons(i),df);
        leadMotorDiff(i).EverythingElse = (temp(1,2)+temp(1,3)) - (temp(1,4)+temp(1,5));
        leadDiffAxons(i).EverythingElse = uniqueEverythingElseAxons(i);
        plot([1,2],[leadDiff(i).EverythingElse,leadMotorSynapseDiff(i).EverythingElse],'-','Color',[1,0,0,0.1]);
        hold on;
    else
        subplot(5,4,17)
        lagDiff(i).EverythingElse = sum(leadEverythingElseAxons == uniqueEverythingElseAxons(i)) - sum(lagEverythingElseAxons == uniqueEverythingElseAxons(i));
        l = find(everythginElseMotorDist(:,1) == uniqueEverythingElseAxons(i));
        lagMotorSynapseDiff(i).EverythingElse = (everythginElseMotorDist(l(1),2)+everythginElseMotorDist(l(1),3)) - (everythginElseMotorDist(l(1),4)+everythginElseMotorDist(l(1),4));
        temp  = isPostSynapseMotor(uniqueEverythingElseAxons(i),df);
        lagMotorDiff(i).EverythingElse = (temp(1,2)+temp(1,3)) - (temp(1,4)+temp(1,5));
        lagDiffAxons(i).EverythingElse = uniqueEverythingElseAxons(i)
        plot([1,2],[lagDiff(i).EverythingElse,lagMotorSynapseDiff(i).EverythingElse],'-','Color',[0,0,1,0.1]);
        hold on;
    end
end

plot([1,2],[mean([leadDiff.EverythingElse]), mean([leadMotorSynapseDiff.EverythingElse])],'-o','Color',[1,0,0]);
plot([1,2],[mean([lagDiff.EverythingElse]), mean([lagMotorSynapseDiff.EverythingElse])],'-o','Color',[0,0,1]);
set(gca,'XLim',[0,3]);
title('Remaining axons');
box off;
axis square;

subplot(5,4,18)
scatter([leadDiff.EverythingElse],[leadMotorSynapseDiff.EverythingElse],20,'o','MarkerFaceColor',colors(5,:),'MarkerEdgeColor','none');
hold on;
scatter([lagDiff.EverythingElse],[lagMotorSynapseDiff.EverythingElse],20,'o','MarkerFaceColor',colors(5,:),'MarkerEdgeColor','none');
line([0,0],[-100,100],'linestyle',':','color','k');
line([-10,10],[0,0],'linestyle',':','color','k');
box off;
axis square;

subplot(5,4,19)

scatter([leadDiff.EverythingElse],[leadMotorDiff.EverythingElse],20,'o','MarkerFaceColor',colors(5,:),'MarkerEdgeColor','none');
hold on;
scatter([lagDiff.EverythingElse],[lagMotorDiff.EverythingElse],20,'o','MarkerFaceColor',colors(5,:),'MarkerEdgeColor','none');
line([0,0],[-20,20],'linestyle',':','color','k');
line([-10,10],[0,0],'linestyle',':','color','k');
ylabel('#(ABD - ABDi) neurons');

box off;
axis square;

subplot(5,4,20)
scatter([leadMotorSynapseDiff.EverythingElse],[leadMotorDiff.EverythingElse],20,'o','MarkerFaceColor',leadColor,'MarkerEdgeColor','none');
hold on
scatter([lagMotorSynapseDiff.EverythingElse],[lagMotorDiff.EverythingElse],20,'o','MarkerFaceColor',lagColor,'MarkerEdgeColor','none');
line([0,0],[-20,20],'linestyle',':','color','k');
line([-200,200],[0,0],'linestyle',':','color','k');
ylabel('#(ABD - ABDi) neurons');
xlabel('#(ABD - ABDi) synapses');
text(150,-15,sprintf('n=%d',size(uniqueEverythingElseAxons,1)));

box off;
axis square;

%save variables

save('leadDiff.mat','leadDiff')
save('leadMotorDiff.mat','leadMotorDiff');
save('leadDiffAxons.mat','leadDiffAxons');

save('lagDiff.mat','lagDiff')
save('lagMotorDiff.mat','lagMotorDiff');
save('lagDiffAxons.mat','lagDiffAxons');

%% Pie plots of input/putput fractions


leadPie = [sum(vertcat(Lead.isSaccadic)),sum(vertcat(Lead.isVestibular)),sum(vertcat(Lead.isIntegrator)), ...
    sum(vertcat(Lead.isContra)),sum(vertcat(Lead.isEverythingElse))];

lagPie = [sum(vertcat(Lag.isSaccadic)),sum(vertcat(Lag.isVestibular)),sum(vertcat(Lag.isIntegrator)), ...
    sum(vertcat(Lag.isContra)),sum(vertcat(Lag.isEverythingElse))];

figure;

subplot(3,5,1)
p1 = pie(leadPie);
txt = {'Saccadic';'Vestibular';'Integrator'; 'Contra';'Rest'}; 
legend(txt,'Location','bestoutside');
colormap(colors(1:5,:));
title('Lead makeup');

subplot(3,5,4)
p2 = pie(lagPie);
txt = {'Saccadic';'Vestibular';'Integrator'; 'Contra';'Rest'}; 
colormap(colors(1:5,:));
title('lag makeup');

%  plot spatial locaiton onto integrator neurons

for i = 1:numel(leadNeurons)    
    leadSaccadicPathlength(1:size(Lead(i).Saccadic,1),i) = Lead(i).PathLength(Lead(i).isSaccadic)/max(Pvec_tree(Lead(i).Tree{1}));
    leadVestibularPathlength(1:size(Lead(i).Vestibular,1),i) = Lead(i).PathLength(Lead(i).isVestibular)/max(Pvec_tree(Lead(i).Tree{1}));
    leadIntegratorPathlength(1:size(Lead(i).Integrator,1),i) = Lead(i).PathLength(Lead(i).isIntegrator)/max(Pvec_tree(Lead(i).Tree{1}));
    leadContraPathlength(1:size(Lead(i).Contra,1),i) = Lead(i).PathLength(Lead(i).isContra)/max(Pvec_tree(Lead(i).Tree{1}));
    leadEverythingElsePathlength(1:size(Lead(i).EverythingElse,1),i) = Lead(i).PathLength(Lead(i).isEverythingElse)/max(Pvec_tree(Lead(i).Tree{1}));
end

% lag neurons

for i = 1:numel(lagNeurons)    
    lagSaccadicPathlength(1:size(Lag(i).Saccadic,1),i) = Lag(i).PathLength(Lag(i).isSaccadic)/max(Pvec_tree(Lag(i).Tree{1}));
    lagVestibularPathlength(1:size(Lag(i).Vestibular,1),i) = Lag(i).PathLength(Lag(i).isVestibular)/max(Pvec_tree(Lag(i).Tree{1}));
    lagIntegratorPathlength(1:size(Lag(i).Integrator,1),i) = Lag(i).PathLength(Lag(i).isIntegrator)/max(Pvec_tree(Lag(i).Tree{1}));
    lagContraPathlength(1:size(Lag(i).Contra,1),i) = Lag(i).PathLength(Lag(i).isContra)/max(Pvec_tree(Lag(i).Tree{1}));
    lagEverythingElsePathlength(1:size(Lag(i).EverythingElse,1),i) = Lag(i).PathLength(Lag(i).isEverythingElse)/max(Pvec_tree(Lag(i).Tree{1}));
end

subplot(3,5,6)
histogram(leadSaccadicPathlength(leadSaccadicPathlength~=0),20,'FaceColor',leadColor,'EdgeColor','none',...
    'Normalization','probability');
hold on;
histogram(lagSaccadicPathlength(lagSaccadicPathlength~=0),20,'FaceColor',lagColor,'EdgeColor','none',...
    'Normalization','probability');
title('Saccadic axons');
ylabel('Probability');
box off;
axis square;

subplot(3,5,7)
histogram(leadVestibularPathlength(leadVestibularPathlength~=0),20,'FaceColor',leadColor,'EdgeColor','none',...
    'Normalization','probability');
hold on;
histogram(lagVestibularPathlength(lagVestibularPathlength~=0),20,'FaceColor',lagColor,'EdgeColor','none',...
    'Normalization','probability');
title('Vestibular axons');

box off;
axis square;

subplot(3,5,8)
histogram(leadIntegratorPathlength(leadIntegratorPathlength~=0),20,'FaceColor',leadColor,'EdgeColor','none',...
    'Normalization','probability');
hold on;
histogram(lagIntegratorPathlength(lagIntegratorPathlength~=0),20,'FaceColor',lagColor,'EdgeColor','none',...
    'Normalization','probability');
title('Integrator axons');

box off;
axis square;

subplot(3,5,9)
histogram(leadContraPathlength(leadContraPathlength~=0),20,'FaceColor',leadColor,'EdgeColor','none',...
    'Normalization','probability');
hold on;
histogram(lagContraPathlength(lagContraPathlength~=0),20,'FaceColor',lagColor,'EdgeColor','none',...
    'Normalization','probability');
title('Contra axons');

box off;
axis square;

subplot(3,5,10)
histogram(leadEverythingElsePathlength(leadEverythingElsePathlength~=0),20,'FaceColor',leadColor,'EdgeColor','none',...
    'Normalization','probability');
hold on;
histogram(lagEverythingElsePathlength(lagEverythingElsePathlength~=0),20,'FaceColor',lagColor,'EdgeColor','none',...
    'Normalization','probability');
title('Rest axons');
xlabel('Pathlength (um)')
box off;
axis square;

% distributions based on rhombomere organization

%r3 organization

for i = 1:numel(leadNeurons)    
    leadR3Pathlength(1:size(Lead(i).r3,1),i) = Lead(i).PathLength(Lead(i).r3(:,2))/max(Pvec_tree(Lead(i).Tree{1}));
    leadR4Pathlength(1:size(Lead(i).r4,1),i) = Lead(i).PathLength(Lead(i).r4(:,2))/max(Pvec_tree(Lead(i).Tree{1}));
    leadR5Pathlength(1:size(Lead(i).r5,1),i) = Lead(i).PathLength(Lead(i).r5(:,2))/max(Pvec_tree(Lead(i).Tree{1}));
    leadR6Pathlength(1:size(Lead(i).r6,1),i) = Lead(i).PathLength(Lead(i).r6(:,2))/max(Pvec_tree(Lead(i).Tree{1}));
    leadR7Pathlength(1:size(Lead(i).r7,1),i) = Lead(i).PathLength(Lead(i).r7(:,2))/max(Pvec_tree(Lead(i).Tree{1}));
end

for i = 1:numel(lagNeurons)  
    lagR3Pathlength(1:size(Lag(i).r3,1),i) = Lag(i).PathLength(Lag(i).r3(:,2))/max(Pvec_tree(Lag(i).Tree{1}));
    lagR4Pathlength(1:size(Lag(i).r4,1),i) = Lag(i).PathLength(Lag(i).r4(:,2))/max(Pvec_tree(Lag(i).Tree{1}));
    lagR5Pathlength(1:size(Lag(i).r5,1),i) = Lag(i).PathLength(Lag(i).r5(:,2))/max(Pvec_tree(Lag(i).Tree{1}));
    lagR6Pathlength(1:size(Lag(i).r6,1),i) = Lag(i).PathLength(Lag(i).r6(:,2))/max(Pvec_tree(Lag(i).Tree{1}));
    lagR7Pathlength(1:size(Lag(i).r7,1),i) = Lag(i).PathLength(Lag(i).r7(:,2))/max(Pvec_tree(Lag(i).Tree{1}));
end

% plot histograms

% subplot(3,5,11)
% histogram(leadR3Pathlength(leadR3Pathlength~=0),20,'FaceColor',leadColor,'EdgeColor','none',...
%     'Normalization','probability');
% hold on;
% histogram(lagR3Pathlength(lagR3Pathlength~=0),20,'FaceColor',lagColor,'EdgeColor','none',...
%     'Normalization','probability');
% title('R3 axons');
% xlabel('Pathlength (um)')
% box off;
% axis square;

subplot(3,5,11)
histogram(leadR3Pathlength(leadR3Pathlength~=0),20,'FaceColor',leadColor,'EdgeColor','none',...
    'Normalization','probability');
hold on;
histogram(lagR3Pathlength(lagR3Pathlength~=0),20,'FaceColor',lagColor,'EdgeColor','none',...
    'Normalization','probability');
title('R3 axons');
xlabel('Pathlength (um)')
box off;
axis square;

subplot(3,5,12)
histogram(leadR4Pathlength(leadR4Pathlength~=0),20,'FaceColor',leadColor,'EdgeColor','none',...
    'Normalization','probability');
hold on;
histogram(lagR4Pathlength(lagR4Pathlength~=0),20,'FaceColor',lagColor,'EdgeColor','none',...
    'Normalization','probability');
title('R4 axons');
xlabel('Pathlength (um)')
box off;
axis square;

subplot(3,5,13)
histogram(leadR5Pathlength(leadR5Pathlength~=0),20,'FaceColor',leadColor,'EdgeColor','none',...
    'Normalization','probability');
hold on;
histogram(lagR5Pathlength(lagR5Pathlength~=0),20,'FaceColor',lagColor,'EdgeColor','none',...
    'Normalization','probability');
title('R5 axons');
xlabel('Pathlength (um)')
box off;
axis square;

subplot(3,5,14)
histogram(leadR6Pathlength(leadR6Pathlength~=0),20,'FaceColor',leadColor,'EdgeColor','none',...
    'Normalization','probability');
hold on;
histogram(lagR6Pathlength(lagR6Pathlength~=0),20,'FaceColor',lagColor,'EdgeColor','none',...
    'Normalization','probability');
title('R6 axons');
xlabel('Pathlength (um)')
box off;
axis square;


subplot(3,5,15)
histogram(leadR7Pathlength(leadR7Pathlength~=0),20,'FaceColor',leadColor,'EdgeColor','none',...
    'Normalization','probability');
hold on;
histogram(lagR7Pathlength(lagR7Pathlength~=0),20,'FaceColor',lagColor,'EdgeColor','none',...
    'Normalization','probability');
title('R7 axons');
xlabel('Pathlength (um)')
box off;
axis square;

%% Looking at unique populations only



% for i = 1:numel(leadNeurons)  
%     var1 = ismember(Lead(i).Inputs,leadOnlySaccadic);
%     var2 = ismember(Lead(i).Inputs,commonSaccadic);
%     leadOnlySaccadicPathlength(1:sum(var1),i) = Lead(i).PathLength(var1)/max(Pvec_tree(Lead(i).Tree{1}));
%     CommonSaccadicPathlength(1:sum(var2),i) = Lead(i).PathLength(var2)/max(Pvec_tree(Lead(i).Tree{1}));
% end
% 
% for i = 1:numel(lagNeurons)  
%     var1 = ismember(Lag(i).Inputs,lagOnlySaccadic);
%     var2 = ismember(Lag(i).Inputs,commonSaccadic);
%     lagOnlySaccadicPathlength(1:sum(var1),i) = Lag(i).PathLength(var1)/max(Pvec_tree(Lag(i).Tree{1}));
%     CommonSaccadicPathlength(1:sum(var2),numel(leadNeurons)+i) = Lag(i).PathLength(var2)/max(Pvec_tree(Lag(i).Tree{1}));
% end
% 
% commonSaccadicPathlength(:,1:5) = leadCommonSaccadicPathlength;
% commonSaccadicPathlength(:,5:9) =  lagCommonSaccadicPathlength;
% 
% histogram(leadOnlySaccadicPathlength(leadOnlySaccadicPathlength~=0),10);
% hold on;
% histogram(lagOnlySaccadicPathlength(leadOnlySaccadicPathlength~=0),10);
% histogram(CommonSaccadicPathlength(CommonSaccadicPathlength~=0),10);
% 



% %% Plot all Lead/Lag axons from connectome
% 
% load ConnMatrixPre.mat;
% load AllCells.mat;
% 
% 
% [~,leadA] = intersect(AllCells,leadDiffAxons);
% [~,leadB] = intersect(AllCells,leadNeurons);
% 
% connMatLeadAxons = zeros(size(leadB,1),size(leadA,1));
% 
% for i = 1:size(leadB,1)
%     connMatLeadAxons(i,:) = ConnMatrixPre(leadB(i),leadA);
% end
% 
% subplot(2,1,1);
% cspy(connMatLeadAxons,'ColorMap',CustomCMap,'Level',9,'MarkerSize',20);
% 
% [~,lagA] = intersect(AllCells, lagDiffAxons);
% [~,lagB] = intersect(AllCells,lagNeurons);
% 
% connMatLagAxons = zeros(size(lagB,1),size(lagA,1));
% 
% for i = 1:size(lagB,1)
%     connMatLagAxons(i,:) = ConnMatrixPre(lagB(i),lagA);
% end
% 
% subplot(2,1,2);
% cspy(connMatLagAxons,'ColorMap',CustomCMap,'Level',9,'MarkerSize',20);
% 
% %% plot the axons on the Z-brain atlas
% 
% figure;
% 
% subplot(1,3,1)
% transform_swc_AV([leadDiffAxons.Saccadic]',[1,0,0],[],false);
% transform_swc_AV([lagDiffAxons.Saccadic]',[0,0,1],[],false);
% title('Saccadic');
% 
% subplot(1,3,2)
% transform_swc_AV([leadDiffAxons.Contra]',[1,0,0],[],false);
% transform_swc_AV([lagDiffAxons.Contra]',[0,0,1],[],false);
% title('Contra.');
% 
% subplot(1,3,3)
% transform_swc_AV([leadDiffAxons.Integrator]',[1,0,0],[],false);
% transform_swc_AV([lagDiffAxons.Integrator]',[0,0,1],[],false);
% title('Integrator');
% 
% 
% %% plot location on cells based on identity
% 
% figure;
% for i = 1:numel(leadNeurons)
% subplot(4,4,i)
% plot_tree(Lead(i).Tree{1},'k',[],[],[],'-3l');
% hold on;
% scatter3(Lead(i).PreSynCoordsTransformed(Lead(i).isSaccadic,1), Lead(i).PreSynCoordsTransformed(Lead(i).isSaccadic,2), Lead(i).PreSynCoordsTransformed(Lead(i).isSaccadic,3),'o',...
%     'MarkerFaceColor',colors(1,:),'MarkerEdgeColor','none')
% scatter3(Lead(i).PreSynCoordsTransformed(Lead(i).isContra,1), Lead(i).PreSynCoordsTransformed(Lead(i).isContra,2), Lead(i).PreSynCoordsTransformed(Lead(i).isContra,3),'o', ...
%      'MarkerFaceColor',colors(2,:),'MarkerEdgeColor','none')
% scatter3(Lead(i).PreSynCoordsTransformed(Lead(i).isIntegrator,1), Lead(i).PreSynCoordsTransformed(Lead(i).isIntegrator,2), Lead(i).PreSynCoordsTransformed(Lead(i).isIntegrator,3),'o', ...
%  'MarkerFaceColor',colors(3,:),'MarkerEdgeColor','none')
% title(leadNeurons(i));
% axis off;
% end
% sgtitle('Lead Neurons');
% 
% for i = 1:numel(lagNeurons)
% subplot(4,4,8+i)
% plot_tree(Lag(i).Tree{1},'k',[],[],[],'-3l');
% hold on;
% scatter3(Lag(i).PreSynCoordsTransformed(Lag(i).isSaccadic,1), Lag(i).PreSynCoordsTransformed(Lag(i).isSaccadic,2), Lag(i).PreSynCoordsTransformed(Lag(i).isSaccadic,3),'o',...
%     'MarkerFaceColor',colors(1,:),'MarkerEdgeColor','none')
% scatter3(Lag(i).PreSynCoordsTransformed(Lag(i).isContra,1), Lag(i).PreSynCoordsTransformed(Lag(i).isContra,2), Lag(i).PreSynCoordsTransformed(Lag(i).isContra,3),'o', ...
%      'MarkerFaceColor',colors(2,:),'MarkerEdgeColor','none')
% scatter3(Lag(i).PreSynCoordsTransformed(Lag(i).isIntegrator,1), Lag(i).PreSynCoordsTransformed(Lag(i).isIntegrator,2), Lag(i).PreSynCoordsTransformed(Lag(i).isIntegrator,3),'o', ...
%  'MarkerFaceColor',colors(3,:),'MarkerEdgeColor','none')
% title(lagNeurons(i));
% axis off;
% end
% sgtitle('Lag Neurons');
% 
% 
% 
% %% plot location on cells based on rhombomere
% 
% figure;
% for i = 1:numel(leadNeurons)
% subplot(4,4,i)
% plot_tree(Lead(i).Tree{1},'k',[],[],[],'-3l');
% hold on;
% %r3
% scatter3(Lead(i).PreSynCoordsTransformed(Lead(i).r3(:,2),1),Lead(i).PreSynCoordsTransformed(Lead(i).r3(:,2),2), ...
%     Lead(i).PreSynCoordsTransformed(Lead(i).r3(:,2),3),'MarkerFaceColor',colors(1,:),'MarkerEdgeColor','none');
% %r4
% scatter3(Lead(i).PreSynCoordsTransformed(Lead(i).r4(:,2),1),Lead(i).PreSynCoordsTransformed(Lead(i).r4(:,2),2), ...
%     Lead(i).PreSynCoordsTransformed(Lead(i).r4(:,2),3),'MarkerFaceColor',colors(2,:),'MarkerEdgeColor','none');
% %r5
% scatter3(Lead(i).PreSynCoordsTransformed(Lead(i).r5(:,2),1),Lead(i).PreSynCoordsTransformed(Lead(i).r5(:,2),2), ...
%     Lead(i).PreSynCoordsTransformed(Lead(i).r5(:,2),3),'MarkerFaceColor',colors(3,:),'MarkerEdgeColor','none');
% %r6
% scatter3(Lead(i).PreSynCoordsTransformed(Lead(i).r6(:,2),1),Lead(i).PreSynCoordsTransformed(Lead(i).r6(:,2),2), ...
%     Lead(i).PreSynCoordsTransformed(Lead(i).r6(:,2),3),'MarkerFaceColor',colors(4,:),'MarkerEdgeColor','none');
% %r7
% scatter3(Lead(i).PreSynCoordsTransformed(Lead(i).r7(:,2),1),Lead(i).PreSynCoordsTransformed(Lead(i).r7(:,2),2), ...
%     Lead(i).PreSynCoordsTransformed(Lead(i).r7(:,2),3),'MarkerFaceColor',colors(5,:),'MarkerEdgeColor','none');
% 
% title(leadNeurons(i));
% axis off;
% end
% 
% sgtitle('Lead Neurons');
% 
% for i = 1:numel(lagNeurons)
% subplot(4,4,8+i)
% plot_tree(Lag(i).Tree{1},'k',[],[],[],'-3l');
% hold on;
% %r3
% scatter3(Lag(i).PreSynCoordsTransformed(Lag(i).r3(:,2),1),Lag(i).PreSynCoordsTransformed(Lag(i).r3(:,2),2), ...
%     Lag(i).PreSynCoordsTransformed(Lag(i).r3(:,2),3),'MarkerFaceColor',colors(1,:),'MarkerEdgeColor','none');
% %r4
% scatter3(Lag(i).PreSynCoordsTransformed(Lag(i).r4(:,2),1),Lag(i).PreSynCoordsTransformed(Lag(i).r4(:,2),2), ...
%     Lag(i).PreSynCoordsTransformed(Lag(i).r4(:,2),3),'MarkerFaceColor',colors(2,:),'MarkerEdgeColor','none');
% %r5
% scatter3(Lag(i).PreSynCoordsTransformed(Lag(i).r5(:,2),1),Lag(i).PreSynCoordsTransformed(Lag(i).r5(:,2),2), ...
%     Lag(i).PreSynCoordsTransformed(Lag(i).r5(:,2),3),'MarkerFaceColor',colors(3,:),'MarkerEdgeColor','none');
% %r6
% scatter3(Lag(i).PreSynCoordsTransformed(Lag(i).r6(:,2),1),Lag(i).PreSynCoordsTransformed(Lag(i).r6(:,2),2), ...
%     Lag(i).PreSynCoordsTransformed(Lag(i).r6(:,2),3),'MarkerFaceColor',colors(4,:),'MarkerEdgeColor','none');
% %r7
% scatter3(Lag(i).PreSynCoordsTransformed(Lag(i).r7(:,2),1),Lag(i).PreSynCoordsTransformed(Lag(i).r7(:,2),2), ...
%     Lag(i).PreSynCoordsTransformed(Lag(i).r7(:,2),3),'MarkerFaceColor',colors(5,:),'MarkerEdgeColor','none');
% 
% title(lagNeurons(i));
% axis off;
% end
% sgtitle('Lag Neurons');
%     
% 
% %% plot based on rhombomere pop
% 
% 
%  %r3
%  leadR3 = vertcat(Lead.r3);
%  lagR3 = vertcat(Lag.r3);
% 
%  transform_swc_AV(leadR3(:,1),[1,0,0],[],false);
%  transform_swc_AV(lagR3(:,1),[0,0,1],[],false);
%  
%   %r4
%  leadR4 = vertcat(Lead.r4);
%  lagR4 = vertcat(Lag.r4);
% 
%  transform_swc_AV(leadR4(:,1),[1,0,0],[],false);
%  transform_swc_AV(lagR4(:,1),[0,0,1],[],false);
%  
%   %r5
%  leadR5 = vertcat(Lead.r5);
%  lagR5 = vertcat(Lag.r5);
% 
%  transform_swc_AV(leadR5(:,1),[1,0,0],[],false);
%  transform_swc_AV(lagR5(:,1),[0,0,1],[],false);
%  
%   %r3
%  leadR6 = vertcat(Lead.r6);
%  lagR6 = vertcat(Lag.r6);
% 
%  transform_swc_AV(leadR6(:,1),[1,0,0],[],false);
%  transform_swc_AV(lagR6(:,1),[0,0,1],[],false);
%  
%   %r3
%  leadR7 = vertcat(Lead.r7);
%  lagR7 = vertcat(Lag.r7);
% 
%  transform_swc_AV(leadR7(:,1),[1,0,0],[],false);
%  transform_swc_AV(lagR7(:,1),[0,0,1],[],false);
% 
 %% make pretty plots
 
 colors2 = colorcet('CBTD1','N',5,'reverse',1); % blue-->red for the partners
 colors1 = colorcet('D1A','N',5, 'reverse',1); % blue -->red for the neurons
 
%colors = colorcet('D1A','N',20,'reverse',1);
%colors = distinguishable_colors(5,{'w','k'});
 
leadOnlySaccadic = setdiff(leadSaccadicAxons, lagSaccadicAxons);
lagOnlySaccadic = setdiff(lagSaccadicAxons, leadSaccadicAxons);
commonSaccadic = intersect(leadSaccadicAxons,lagSaccadicAxons);
 
 for i = 1:numel(commonSaccadic)
     commonSaccadicDiff(i).ID = commonSaccadic(i);
     commonSaccadicDiff(i).Lead = sum(ismember(vertcat(Lead.Saccadic), commonSaccadic(i)));
     commonSaccadicDiff(i).Lag = sum(ismember(vertcat(Lag.Saccadic), commonSaccadic(i)));
 end
 
 index = find([commonSaccadicDiff.Lead]>[commonSaccadicDiff.Lag]);
 leadOnlySaccadic = [leadOnlySaccadic;[commonSaccadicDiff(index).ID]'];
 index = find([commonSaccadicDiff.Lead]<[commonSaccadicDiff.Lag]);
 lagOnlySaccadic = [lagOnlySaccadic;[commonSaccadicDiff(index).ID]'];

 % plot lead and Lag neurons
 subplot(1,2,1)
 transform_swc_AV(leadNeurons,colors1(1:2,:),[],true,false);
 transform_swc_AV(leadOnlySaccadic,colors2(1:2,:),[],false,false);
 
%export_fig('/Users/ashwin/Desktop/SaccadicsLEADOntoInt.png','-r300','-transparent');
%close all;

 subplot(1,2,2)
 transform_swc_AV(lagNeurons,colors1(4:5,:),[],true,false);
 transform_swc_AV(lagOnlySaccadic,colors2(4:5,:),[],false,false);
 

 export_fig('/Users/ashwin/Desktop/SaccadicsLAGOntoInt.png','-r300','-transparent');
 close all;
% 
%  %transform_swc_AV(commonSaccadic,[colors(9,:)+colors(12,:)]/3,[],false,false);
% 
% % vestibular pop
% 
%  
% leadOnlyVestibular = setdiff(leadVestibularAxons, lagVestibularAxons);
% lagOnlyVestibular  = setdiff(lagVestibularAxons, leadVestibularAxons);
% commonVestibular  = intersect(leadVestibularAxons,lagVestibularAxons);
% 
% for i = 1:numel(commonVestibular)
%      commonVestibularDiff(i).ID = commonVestibular(i);
%      commonVestibularDiff(i).Lead = sum(ismember(vertcat(Lead.Vestibular), commonVestibular(i)));
%      commonVestibularDiff(i).Lag = sum(ismember(vertcat(Lag.Vestibular), commonVestibular(i)));
%  end
%  
%  index = find([commonVestibularDiff.Lead]>[commonVestibularDiff.Lag]);
%  leadOnlyVestibular = [leadOnlyVestibular;[commonVestibularDiff(index).ID]'];
%  index = find([commonVestibularDiff.Lead]<[commonVestibularDiff.Lag]);
%  lagOnlyVestibular = [lagOnlyVestibular;[commonVestibularDiff(index).ID]'];
% 
%  % plot lead and Lag neurons
%  subplot(1,2,1)
%  transform_swc_AV(leadNeurons,colors(9,:),[],true,false);
%  transform_swc_AV(leadOnlyVestibular,colors(4,:),[],false,false);
% export_fig('/Users/ashwin/Desktop/VestibularLEADOntoInt.png','-r300','-transparent');
% close all;
%  
%  subplot(1,2,2)
%  transform_swc_AV(lagNeurons,colors(12,:),[],true,false);
%  transform_swc_AV(lagOnlyVestibular,colors(16,:),[],false,false);
%  export_fig('/Users/ashwin/Desktop/VestibularLAGOntoInt.png','-r300','-transparent');
% close all;
%  
%  % Contra pop
% 
%  
% leadOnlyContra = setdiff(leadContraAxons, lagContraAxons);
% lagOnlyContra  = setdiff(lagContraAxons, leadContraAxons);
% commonContra  = intersect(leadContraAxons,lagContraAxons);
% 
% for i = 1:numel(commonContra)
%      commonContraDiff(i).ID = commonContra(i);
%      commonContraDiff(i).Lead = sum(ismember(vertcat(Lead.Contra), commonContra(i)));
%      commonContraDiff(i).Lag = sum(ismember(vertcat(Lag.Contra), commonContra(i)));
%  end
%  
%  index = find([commonContraDiff.Lead]>[commonContraDiff.Lag]);
%  leadOnlyContra = [leadOnlyContra;[commonContraDiff(index).ID]'];
%  index = find([commonContraDiff.Lead]<[commonContraDiff.Lag]);
%  lagOnlyContra = [lagOnlyContra;[commonContraDiff(index).ID]'];
% 
%  % plot lead and Lag neurons
%  subplot(1,2,1)
%  transform_swc_AV(leadNeurons,colors(9,:),[],true,false);
%  transform_swc_AV(leadOnlyContra,colors(4,:),[],false,false);
%  export_fig('/Users/ashwin/Desktop/ContraLEADOntoInt.png','-r300','-transparent');
% close all;
%  
%  subplot(1,2,2)
%  transform_swc_AV(lagNeurons,colors(12,:),[],true,false);
%  transform_swc_AV(lagOnlyContra,colors(16,:),[],false,false);
%   export_fig('/Users/ashwin/Desktop/ContraLAGOntoInt.png','-r300','-transparent');
% close all;
%  
% % seperate the Contra by R5/R6 and remaining
% clear temp;
% R4LeadContra = [];
% R5LeadContra = [];
% R6LeadContra = [];
% 
% for i = 1:numel(leadNeurons)
%     temp = find(ismember(Lead(i).Inputs,leadOnlyContra));
%     R4LeadContra = [R4LeadContra;Lead(i).InputsRhombomeres(temp(find(Lead(i).InputsRhombomeres(temp,6) ==1)),1)];
%     R5LeadContra = [R5LeadContra;Lead(i).InputsRhombomeres(temp(find(Lead(i).InputsRhombomeres(temp,7) ==1)),1)];
%     R6LeadContra = [R6LeadContra;Lead(i).InputsRhombomeres(temp(find(Lead(i).InputsRhombomeres(temp,8) ==1)),1)];
%     clear temp;
% end
% 
% 
% 
% 
% 
% clear temp;
% R4LagContra = [];
% R5LagContra = [];
% R6LagContra = [];
% 
% for i = 1:numel(lagNeurons)
%     temp = find(ismember(Lag(i).Inputs,lagOnlyContra));
%     R4LagContra = [R4LagContra;Lag(i).InputsRhombomeres(temp(find(Lag(i).InputsRhombomeres(temp,6) ==1)),1)];
%     R5LagContra = [R5LagContra;Lag(i).InputsRhombomeres(temp(find(Lag(i).InputsRhombomeres(temp,7) ==1)),1)];
%     R6LagContra = [R6LagContra;Lag(i).InputsRhombomeres(temp(find(Lag(i).InputsRhombomeres(temp,8) ==1)),1)];
%     clear temp;
% end
% 
% subplot(1,3,1)
% transform_swc_AV(R4LeadContra,colors(4,:),[],true,false);
% transform_swc_AV(R4LagContra,colors(16,:),[],false,false);
% 
% subplot(1,3,2)
% transform_swc_AV(R5LeadContra,colors(4,:),[],true,false);
% transform_swc_AV(R5LagContra,colors(16,:),[],false,false);
% 
% subplot(1,3,3)
% transform_swc_AV(R6LeadContra,colors(4,:),[],true,false);
% transform_swc_AV(R6LagContra,colors(16,:),[],false,false);




