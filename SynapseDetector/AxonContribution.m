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
    leadSaccadicPathlength(1:size(Lead(i).Saccadic,1),i) = Lead(i).PathLength(Lead(i).isSaccadic);
    leadVestibularPathlength(1:size(Lead(i).Vestibular,1),i) = Lead(i).PathLength(Lead(i).isVestibular);
    leadIntegratorPathlength(1:size(Lead(i).Integrator,1),i) = Lead(i).PathLength(Lead(i).isIntegrator);
    leadContraPathlength(1:size(Lead(i).Contra,1),i) = Lead(i).PathLength(Lead(i).isContra);
    leadEverythingElsePathlength(1:size(Lead(i).EverythingElse,1),i) = Lead(i).PathLength(Lead(i).isEverythingElse);
end

% lag neurons

for i = 1:numel(lagNeurons)    
    lagSaccadicPathlength(1:size(Lag(i).Saccadic,1),i) = Lag(i).PathLength(Lag(i).isSaccadic);
    lagVestibularPathlength(1:size(Lag(i).Vestibular,1),i) = Lag(i).PathLength(Lag(i).isVestibular);
    lagIntegratorPathlength(1:size(Lag(i).Integrator,1),i) = Lag(i).PathLength(Lag(i).isIntegrator);
    lagContraPathlength(1:size(Lag(i).Contra,1),i) = Lag(i).PathLength(Lag(i).isContra);
    lagEverythingElsePathlength(1:size(Lag(i).EverythingElse,1),i) = Lag(i).PathLength(Lag(i).isEverythingElse);
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




%%  plot axon classes onto Motor population
%saccadic
leadSaccadicOnMotor = zeros(size(all.Saccadic,1),5);
lagSaccadicOnMotor = zeros(size(all.Saccadic,1),5);


for i =1:size(allSaccadic,1)
    l = find(allSaccadic(i) == allLeadSaccadicMotorDist(:));
    if l~=0
        leadSaccadicOnMotor(i,:) = allLeadSaccadicMotorDist(l,:);
        clear l;
    end
    
    l = find(allSaccadic(i) == allLagSaccadicMotorDist(:));
    if l~=0
        lagSaccadicOnMotor(i,:) = allLagSaccadicMotorDist(l,:);
        clear l;
    end
end

CustomLeadMap =  cbrewer('seq','Reds',10);
CustomLagMap =  cbrewer('seq','Blues',10);

cmax = max(max([leadSaccadicOnMotor(:,2:end);lagSaccadicOnMotor(:,2:end)]));
cmin = 0;
subplot(4,4,5)
image(leadSaccadicOnMotor(:,2:end),'CDataMapping','scaled');
colormap(gca,CustomLeadMap);
caxis([cmin,cmax]);
colorbar;
xticks([1,2,3,4]);
xticklabels({'ABD_r','ABD_c','ABDI_r','ADBI_c'});
xtickangle(45);
ylabel('Saccadic Axons');
axis square;
box off;
daspect([1,10,4]);

subplot(4,4,6)
image(lagSaccadicOnMotor(:,2:end),'CDataMapping','scaled');
colormap(gca,CustomLagMap);
caxis([cmin,cmax]);
colorbar
axis square;
box off;
daspect([1,10,2]);

% Vestibular
leadVestibularOnMotor = zeros(size(allVest,1),5);
lagVestibularOnMotor = zeros(size(allVest,1),5);

for i =1:size(allVest,1)
    l = find(allVest(i) == allLeadVestibularMotorDist(:));
    if l~=0
        leadVestibularOnMotor(i,:) = allLeadVestibularMotorDist(l,:);
        clear l;
    end
    
    l = find(allVest(i) == allLagVestibularMotorDist(:));
    if l~=0
        lagVestibularOnMotor(i,:) = allLagVestibularMotorDist(l,:);
        clear l;
    end
end

cmax = max(max([leadVestibularOnMotor(:,2:end);lagVestibularOnMotor(:,2:end)]));
cmin = 0;

subplot(4,4,7)
image(leadVestibularOnMotor(:,2:end),'CDataMapping','scaled');
colormap(gca,CustomLeadMap)
colorbar;
xticks([1,2,3,4]);
xticklabels({'ABD_r','ABD_c','ABDI_r','ADBI_c'});
xtickangle(45);
ylabel('Vestibular Axons');
axis square;
box off;
daspect([1,10,4]);

subplot(4,4,8)
image(lagVestibularOnMotor(:,2:end),'CDataMapping','scaled');
colormap(gca,CustomLagMap);
colorbar
axis square;
box off;
daspect([1,10,2]);

% contra
leadContraOnMotor = zeros(size(allSaccadic,1),5);
lagContraOnMotor = zeros(size(allSaccadic,1),5);


for i =1:size(allContra,1)
    l = find(allContra(i) == allLeadContraMotorDist(:));
    if l~=0
        leadContraOnMotor(i,:) = allLeadContraMotorDist(l,:);
        clear l;
    end
    
    l = find(allContra(i) == allLagContraMotorDist(:));
    if l~=0
        lagContraOnMotor(i,:) = allLagContraMotorDist(l,:);
        clear l;
    end
end

cmax = max(max([leadContraOnMotor(:,2:end);lagContraOnMotor(:,2:end)]));
cmin = 0;
subplot(4,4,9)
image(leadContraOnMotor(:,2:end),'CDataMapping','scaled');
colormap(gca,CustomLeadMap);
caxis([cmin,cmax]);
colorbar;
xticks([1,2,3,4]);
xticklabels({'ABD_r','ABD_c','ABDI_r','ADBI_c'});
xtickangle(45);
ylabel('Contra Axons');
axis square;
box off;
daspect([1,10,4]);

subplot(4,4,10)
image(lagContraOnMotor(:,2:end),'CDataMapping','scaled');
colormap(gca,CustomLagMap);
caxis([cmin,cmax]);
colorbar
axis square;
box off;
daspect([1,10,4]);

% integrators
leadIntegratorOnMotor = zeros(size(allIntegrator,1),5);
lagIntegraotorOnMotor = zeros(size(allIntegrator,1),5);

for i =1:size(allIntegrator,1)
    l = find(allIntegrator(i) == allLeadIntegratorMotorDist(:));
    if l~=0
        leadIntegratorOnMotor(i,:) = allLeadIntegratorMotorDist(l,:);
        clear l;
    end
    
    l = find(allIntegrator(i) == allLagIntegratorMotorDist(:));
    if l~=0
        lagIntegraotorOnMotor(i,:) = allLagIntegratorMotorDist(l,:);
        clear l;
    end
end

cmax = max(max([leadIntegratorOnMotor(:,2:end);lagIntegraotorOnMotor(:,2:end)]));
cmin = 0;

subplot(4,4,11)
image(leadIntegratorOnMotor(:,2:end),'CDataMapping','scaled');
colormap(gca,CustomLeadMap)
colorbar;
xticks([1,2,3,4]);
xticklabels({'ABD_r','ABD_c','ABDI_r','ADBI_c'});
xtickangle(45);
ylabel('Integrator Axons');
axis square;
box off;
daspect([1,10,4]);

subplot(4,4,12)
%imagesc(lagEverythingElseOnMotor(:,2:end));
image(lagIntegraotorOnMotor(:,2:end),'CDataMapping','scaled')
colormap(gca,CustomLagMap);
colorbar
axis square;
box off;
daspect([1,10,4]);

% all other axons

leadEverythingElseOnMotor = zeros(size(allRest,1),5);
lagEverythingElseOnMotor = zeros(size(allRest,1),5);

for i =1:size(allRest,1)
    l = find(allRest(i) == allRemainingLeadMotorDist(:));
    if l~=0
        leadEverythingElseOnMotor(i,:) = allRemainingLeadMotorDist(l,:);
        clear l;
    end
    
    l = find(allRest(i) == allRemainingLagMotorDist(:));
    if l~=0
        lagEverythingElseOnMotor(i,:) = allRemainingLagMotorDist(l,:);
        clear l;
    end
end

cmax = max(max([leadEverythingElseOnMotor(:,2:end);lagEverythingElseOnMotor(:,2:end)]));
cmin = 0;

subplot(4,4,13)
image(leadEverythingElseOnMotor(:,2:end),'CDataMapping','scaled');
colormap(gca,CustomLeadMap)
colorbar;
xticks([1,2,3,4]);
xticklabels({'ABD_r','ABD_c','ABDI_r','ADBI_c'});
xtickangle(45);
ylabel('Remaining Axons');
axis square;
box off;
daspect([1,10,4]);

subplot(4,4,14)
%imagesc(lagEverythingElseOnMotor(:,2:end));
image(lagEverythingElseOnMotor(:,2:end),'CDataMapping','scaled')
colormap(gca,CustomLagMap);
colorbar
axis square;
box off;
daspect([1,10,4]);

%%

motorGroups = [1,2,3,4];

subplot(4,4,15);
plot(motorGroups,sum(leadSaccadicOnMotor(:,2:end)),'Color',colors(1,:),'LineWidth',2);
hold on;
plot(motorGroups,sum(leadVestibularOnMotor(:,2:end)),'Color',colors(2,:),'LineWidth',2);
plot(motorGroups,sum(leadContraOnMotor(:,2:end)),'Color',colors(3,:),'LineWidth',2);
plot(motorGroups,sum(leadIntegratorOnMotor(:,2:end)),'Color',colors(4,:),'LineWidth',2);
plot(motorGroups,sum(leadEverythingElseOnMotor(:,2:end)),'Color',colors(5,:),'LineWidth',2);

ylabel('Synapses');
legend({'Sacacde','Vestibular','Contra','Integrator','Rest'},'Location','bestoutside', ...
    'FontSize',7);
axis square;
box off;
set(gca,'XTick',motorGroups,'XTickLabel',{'ABD_r','ABD_c','ABDI_r','ADBI_c'});
set(gca,'color',[1,0,0,0.1]);


subplot(4,4,16)
plot(motorGroups,sum(lagSaccadicOnMotor(:,2:end)),'Color',colors(1,:),'LineWidth',2,'LineStyle','--');
hold on;
plot(motorGroups,sum(lagVestibularOnMotor(:,2:end)),'Color',colors(2,:),'LineWidth',2,'LineStyle','--');
plot(motorGroups,sum(lagContraOnMotor(:,2:end)),'Color',colors(3,:),'LineWidth',2,'LineStyle','--');
plot(motorGroups,sum(lagIntegraotorOnMotor(:,2:end)),'Color',colors(4,:),'LineWidth',2,'LineStyle','--');
plot(motorGroups,sum(lagEverythingElseOnMotor(:,2:end)),'Color',colors(5,:),'LineWidth',2,'LineStyle','--');

ylabel('Synapses');
axis square;
box off;
set(gca,'XTick',motorGroups,'XTickLabel',{'ABD_r','ABD_c','ABDI_r','ADBI_c'},'color',[0,0,1,0.1]);



%% Plot all Lead/Lag axons from connectome

load ConnMatrixPre.mat;
load AllCells.mat;


[~,leadA] = intersect(AllCells, leadDiffAxons);
[~,leadB] = intersect(AllCells,leadNeurons);

connMatLeadAxons = zeros(size(leadB,1),size(leadA,1));

for i = 1:size(leadB,1)
    connMatLeadAxons(i,:) = ConnMatrixPre(leadB(i),leadA);
end

subplot(2,1,1);
cspy(connMatLeadAxons,'ColorMap',CustomCMap,'Level',9,'MarkerSize',20);

[~,lagA] = intersect(AllCells, lagDiffAxons);
[~,lagB] = intersect(AllCells,lagNeurons);

connMatLagAxons = zeros(size(lagB,1),size(lagA,1));

for i = 1:size(lagB,1)
    connMatLagAxons(i,:) = ConnMatrixPre(lagB(i),lagA);
end

subplot(2,1,2);
cspy(connMatLagAxons,'ColorMap',CustomCMap,'Level',9,'MarkerSize',20);

%% plot the axons on the Z-brain atlas
figure(3);

subplot(2,2,1)
transform_swc_AV(allLeadSaccadicMotorDist(:,1)',colors(1,:),colors(1,:),[],false);
subplot(2,2,2);
transform_swc_AV(allLeadVestibularMotorDist(:,1)',colors(2,:),colors(2,:),[],false);
subplot(2,2,3);
transform_swc_AV(allLeadContraMotorDist(:,1)',colors(3,:),colors(3,:),[],false);
subplot(2,2,4);
transform_swc_AV(allRemainingLeadMotorDist(:,1)',colors(4,:),colors(4,:),[],false);




% % plot the location of the synapses
%
% leadSynapseLocations = vertcat(Lead.SynapseLocationTransfromed);
% lagSynapseLocations = vertcat(Lag.SynapseLocationTransfromed);
%
% subplot(1,3,3);
% transform_swc_AV(leadNeurons,colors(5,:),colors(6,:),[],false);
% hold on;
% transform_swc_AV(lagNeurons,colors(1,:),colors(2,:),[],false);
% scatter3(leadSynapseLocations(:,1), leadSynapseLocations(:,2),leadSynapseLocations(:,3),...
%     8,'o','MarkerFaceColor',leadColor,'MarkerEdgeColor','none');
% hold on;
% scatter3(lagSynapseLocations(:,1), lagSynapseLocations(:,2),lagSynapseLocations(:,3),...
%     8,'o','MarkerFaceColor',lagColor,'MarkerEdgeColor','none');
%
%
% daspect([1,1,1]);
% set(gca, 'BoxStyle','full','YDir','reverse','ZDir','reverse');
% set(gca,'ZLim',[0 276],'XLim',[0 495.558], 'YLim',[0 1121.988]); % numbers are from transfrom_swc.m
%

%% Downstream partners of Lead/Lag axons


