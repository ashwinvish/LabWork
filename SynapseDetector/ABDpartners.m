clc;
clear;

if ismac
    addpath(genpath('/Users/ashwin/Documents/LabWork/SynapseDetector/Scripts'));
    df = readtable('/Users/ashwin/Documents/SynapseDetector/11252018.csv');
    fname  = '/Users/ashwin/Documents/LowEMtoHighEM/SWC_all/consensus-20180920/swc/';

else
    addpath(genpath('/usr/people/ashwinv/seungmount/research/Ashwin/Scripts'));
    df = readtable('/usr/people/ashwinv/seungmount/research/Ashwin/SynapseDetector/09202018.csv');
    %ml = readtable('/usr/people/ashwinv/seungmount/research/Ashwin/SynapseDetector/ManualSortList-102218.csv');
    fname =  '/usr/people/ashwinv/seungmount/research/Ashwin/Z-Brain/LowEMtoHighEM/SWC_all/consensus-20180920/swc/'

end

colors = cbrewer('qual','Dark2',10);    

ABDr_CellIDs = [77648, 77710, 77300, 77705, 77305, 77301, 77709, 77672, 77302];
ABDc_CellIDs = [77154, 77646, 77682 ,77628 ,77295 , 77652 ,77292 ,77688 ,77654 ,77658 ,77657 ,77662, 77296];
ABDIr_CellIDs = [77631, 77150, 77618, 77886, 78547, 77158, 78556, 78552, 77665, 77668, 77634, 78553];
ABDIc_CellIDs = [77148, 77625, 77641, 77692, 77144, 77643, 77640, 79051, 79066, 78574];

% complete recontrustucon for 
% ABDr - 77710, 77648
% ABDIr - 78556
% ABDIc - 77641, 77148
%

%%

parfor i = 1:numel(ABDr_CellIDs)
    ABDr(i) = InputsByClass(ABDr_CellIDs(i),df);  
end

parfor i = 1:numel(ABDc_CellIDs)
    ABDc(i) = InputsByClass(ABDc_CellIDs(i),df);
end

parfor i = 1:numel(ABDIr_CellIDs)
    ABDIr(i) = InputsByClass(ABDIr_CellIDs(i),df);
end

parfor i = 1:numel(ABDIc_CellIDs)
    ABDIc(i) = InputsByClass(ABDIc_CellIDs(i),df);
end


%% plot basic statistics

for i = 1:numel(ABDr_CellIDs)
    ABDr(i).NumberofInputs = size(ABDr(i).Inputs,1);
    ABDr(i).NumberOfSaccadic = size(ABDr(i).Saccadic,1);
    ABDr(i).NumberOfVestibular = size(ABDr(i).Vestibular,1);
    ABDr(i).NumberOfContra = size(ABDr(i).Contra,1);
    ABDr(i).NumberOfIntegrator = size(ABDr(i).Integrator,1);
    ABDr(i).NumberOfRemaining = size(ABDr(i).EverythingElse,1);
   
end

for i = 1:numel(ABDc_CellIDs)
    ABDc(i).NumberofInputs = size(ABDc(i).Inputs,1);
    ABDc(i).NumberOfSaccadic = size(ABDc(i).Saccadic,1);
    ABDc(i).NumberOfVestibular = size(ABDc(i).Vestibular,1);
    ABDc(i).NumberOfContra = size(ABDc(i).Contra,1);
    ABDc(i).NumberOfIntegrator = size(ABDc(i).Integrator,1);
    ABDc(i).NumberOfRemaining = size(ABDc(i).EverythingElse,1);
end

for i = 1:numel(ABDIr_CellIDs)
    ABDIr(i).NumberofInputs = size(ABDIr(i).Inputs,1);
    ABDIr(i).NumberOfSaccadic = size(ABDIr(i).Saccadic,1);
    ABDIr(i).NumberOfVestibular = size(ABDIr(i).Vestibular,1);
    ABDIr(i).NumberOfContra = size(ABDIr(i).Contra,1);
    ABDIr(i).NumberOfIntegrator = size(ABDIr(i).Integrator,1);
    ABDIr(i).NumberOfRemaining = size(ABDc(i).EverythingElse,1);
end

for i = 1:numel(ABDIc_CellIDs)
    ABDIc(i).NumberofInputs = size(ABDIc(i).Inputs,1);
    ABDIc(i).NumberOfSaccadic = size(ABDIc(i).Saccadic,1);
    ABDIc(i).NumberOfVestibular = size(ABDIc(i).Vestibular,1);
    ABDIc(i).NumberOfContra = size(ABDIc(i).Contra,1);
    ABDIc(i).NumberOfIntegrator = size(ABDIc(i).Integrator,1);
    ABDIc(i).NumberOfRemaining = size(ABDIc(i).EverythingElse,1);
end

%%
motorGroups = [1,2,3,4];

subplot(4,4,1)
plot(ones(1,size([ABDr.NumberofInputs],2)),[ABDr.NumberofInputs],'ko');
hold on;
plot(2*ones(1,size([ABDc.NumberofInputs],2)),[ABDc.NumberofInputs],'ko');
plot(3*ones(1,size([ABDIr.NumberofInputs],2)),[ABDIr.NumberofInputs],'ko')
plot(4*ones(1,size([ABDIc.NumberofInputs],2)),[ABDIc.NumberofInputs],'ko')

plot([1,2,3,4], [mean([ABDr.NumberofInputs]),mean([ABDc.NumberofInputs]), mean([ABDIr.NumberofInputs]), ...
    mean([ABDIc.NumberofInputs])], '-ko', 'LineWidth',2);

axis square;
box off
set(gca, 'XTick',[1,2,3,4], 'XTickLabel',{'ABDr', 'ABDc', 'ABDIr', 'ABDIc'}, 'XTickLabelRotation',45);
ylabel('total Synapses');

subplot(4,4,2)
plot(ones(1,size([ABDr.NumberOfSaccadic],2)),[ABDr.NumberOfSaccadic],'ko');
hold on;
plot(2*ones(1,size([ABDc.NumberOfSaccadic],2)),[ABDc.NumberOfSaccadic],'ko');
plot(3*ones(1,size([ABDIr.NumberOfSaccadic],2)),[ABDIr.NumberOfSaccadic],'ko')
plot(4*ones(1,size([ABDIc.NumberOfSaccadic],2)),[ABDIc.NumberOfSaccadic],'ko')

plot([1,2,3,4], [mean([ABDr.NumberOfSaccadic]),mean([ABDc.NumberOfSaccadic]), mean([ABDIr.NumberOfSaccadic]), ...
    mean([ABDIc.NumberOfSaccadic])], '-ko', 'LineWidth',2);

axis square;
box off
set(gca, 'XTick',[1,2,3,4], 'XTickLabel',{'ABDr', 'ABDc', 'ABDIr', 'ABDIc'}, 'XTickLabelRotation',45);
ylabel('Saccadic Synapses');


subplot(4,4,3)
plot(ones(1,size([ABDr.NumberOfVestibular],2)),[ABDr.NumberOfVestibular],'ko');
hold on;
plot(2*ones(1,size([ABDc.NumberOfVestibular],2)),[ABDc.NumberOfVestibular],'ko');
plot(3*ones(1,size([ABDIr.NumberOfVestibular],2)),[ABDIr.NumberOfVestibular],'ko')
plot(4*ones(1,size([ABDIc.NumberOfVestibular],2)),[ABDIc.NumberOfVestibular],'ko')

plot([1,2,3,4], [mean([ABDr.NumberOfVestibular]),mean([ABDc.NumberOfVestibular]), mean([ABDIr.NumberOfVestibular]), ...
    mean([ABDIc.NumberOfVestibular])], '-ko', 'LineWidth',2);

axis square;
box off
set(gca, 'XTick',[1,2,3,4], 'XTickLabel',{'ABDr', 'ABDc', 'ABDIr', 'ABDIc'}, 'XTickLabelRotation',45);
ylabel('Vestibular Synapses');


subplot(4,4,4)
plot(ones(1,size([ABDr.NumberOfContra],2)),[ABDr.NumberOfContra],'ko');
hold on;
plot(2*ones(1,size([ABDc.NumberOfContra],2)),[ABDc.NumberOfContra],'ko');
plot(3*ones(1,size([ABDIr.NumberOfContra],2)),[ABDIr.NumberOfContra],'ko')
plot(4*ones(1,size([ABDIc.NumberOfContra],2)),[ABDIc.NumberOfContra],'ko')

plot([1,2,3,4], [mean([ABDr.NumberOfContra]),mean([ABDc.NumberOfContra]), mean([ABDIr.NumberOfContra]), ...
    mean([ABDIc.NumberOfContra])], '-ko', 'LineWidth',2);

axis square;
box off
set(gca, 'XTick',[1,2,3,4], 'XTickLabel',{'ABDr', 'ABDc', 'ABDIr', 'ABDIc'}, 'XTickLabelRotation',45);
ylabel('Contra Synapses');

subplot(4,4,5)
plot(ones(1,size([ABDr.NumberOfIntegrator],2)),[ABDr.NumberOfIntegrator],'ko');
hold on;
plot(2*ones(1,size([ABDc.NumberOfIntegrator],2)),[ABDc.NumberOfIntegrator],'ko');
plot(3*ones(1,size([ABDIr.NumberOfIntegrator],2)),[ABDIr.NumberOfIntegrator],'ko')
plot(4*ones(1,size([ABDIc.NumberOfIntegrator],2)),[ABDIc.NumberOfIntegrator],'ko')

plot([1,2,3,4], [mean([ABDr. NumberOfIntegrator]),mean([ABDc.NumberOfIntegrator]), mean([ABDIr.NumberOfIntegrator]), ...
    mean([ABDIc.NumberOfIntegrator])], '-ko', 'LineWidth',2);

axis square;
box off
set(gca, 'XTick',[1,2,3,4], 'XTickLabel',{'ABDr', 'ABDc', 'ABDIr', 'ABDIc'}, 'XTickLabelRotation',45);
ylabel('Integrator Synapses');

subplot(4,4,6)
plot(ones(1,size([ABDr.NumberOfRemaining],2)),[ABDr.NumberOfRemaining],'ko');
hold on;
plot(2*ones(1,size([ABDc.NumberOfRemaining],2)),[ABDc.NumberOfRemaining],'ko');
plot(3*ones(1,size([ABDIr.NumberOfRemaining],2)),[ABDIr.NumberOfRemaining],'ko')
plot(4*ones(1,size([ABDIc.NumberOfRemaining],2)),[ABDIc.NumberOfRemaining],'ko')

plot([1,2,3,4], [mean([ABDr.NumberOfRemaining]),mean([ABDc.NumberOfRemaining]), mean([ABDIr.NumberOfRemaining]), ...
    mean([ABDIc.NumberOfRemaining])], '-ko', 'LineWidth',2);

axis square;
box off
set(gca, 'XTick',[1,2,3,4], 'XTickLabel',{'ABDr', 'ABDc', 'ABDIr', 'ABDIc'}, 'XTickLabelRotation',45);
ylabel('Remaining Synapses');



%% plot positon of synapses on Z-brain space

% ABDr
for i = 1:numel(ABDr_CellIDs)
    ABDr(i).SaccadicDist = PathLengthToCoordinate(ABDr(i).PreSynCoordsTransformed(ABDr(i).isSaccadic,:),ABDr(i).Tree{1});
    ABDr(i).VestibularDist = PathLengthToCoordinate(ABDr(i).PreSynCoordsTransformed(ABDr(i).isVestibular,:),ABDr(i).Tree{1});
    ABDr(i).ContraDist = PathLengthToCoordinate(ABDr(i).PreSynCoordsTransformed(ABDr(i).isContra,:),ABDr(i).Tree{1}); 
    ABDr(i).IntegratorDist = PathLengthToCoordinate(ABDr(i).PreSynCoordsTransformed(ABDr(i).isIntegrator,:),ABDr(i).Tree{1});
    ABDr(i).EverytingElseDist = PathLengthToCoordinate(ABDr(i).PreSynCoordsTransformed(ABDr(i).isEverythingElse,:), ABDr(i).Tree{1});
end

% ABDc
for i = 1:numel(ABDc_CellIDs)
    ABDc(i).SaccadicDist  = PathLengthToCoordinate(ABDc(i).PreSynCoordsTransformed(ABDc(i).isSaccadic,:),ABDc(i).Tree{1});
    ABDc(i).VestibularDist  = PathLengthToCoordinate(ABDc(i).PreSynCoordsTransformed(ABDc(i).isVestibular,:),ABDc(i).Tree{1});
    ABDc(i).ContraDist  = PathLengthToCoordinate(ABDc(i).PreSynCoordsTransformed(ABDc(i).isContra,:),ABDc(i).Tree{1});
    ABDc(i).IntegratorDist  =  PathLengthToCoordinate(ABDc(i).PreSynCoordsTransformed(ABDc(i).isIntegrator,:),ABDc(i).Tree{1});
    ABDc(i).EverytingElseDist = PathLengthToCoordinate(ABDc(i).PreSynCoordsTransformed(ABDc(i).isEverythingElse,:),ABDc(i).Tree{1});
end

% ABDIr

for i = 1:numel(ABDIr_CellIDs)
    ABDIr(i).SaccadicDist  = PathLengthToCoordinate(ABDIr(i).PreSynCoordsTransformed(ABDIr(i).isSaccadic,:),ABDIr(i).Tree{1});
    ABDIr(i).VestibularDist  = PathLengthToCoordinate(ABDIr(i).PreSynCoordsTransformed(ABDIr(i).isVestibular,:),ABDIr(i).Tree{1});
    ABDIr(i).ContraDist  = PathLengthToCoordinate(ABDIr(i).PreSynCoordsTransformed(ABDIr(i).isContra,:),ABDIr(i).Tree{1});
    ABDIr(i).IntegratorDist  = PathLengthToCoordinate(ABDIr(i).PreSynCoordsTransformed(ABDIr(i).isIntegrator,:),ABDIr(i).Tree{1});
    ABDIr(i).EverytingElseDist = PathLengthToCoordinate(ABDIr(i).PreSynCoordsTransformed(ABDIr(i).isEverythingElse,:),ABDIr(i).Tree{1});
end


% ABDIc

for i = 1:numel(ABDIc_CellIDs)
    ABDIc(i).SaccadicDist  = PathLengthToCoordinate(ABDIc(i).PreSynCoordsTransformed(ABDIc(i).isSaccadic,:),ABDIc(i).Tree{1});
    ABDIc(i).VestibularDist  = PathLengthToCoordinate(ABDIc(i).PreSynCoordsTransformed(ABDIc(i).isVestibular,:),ABDIc(i).Tree{1});
    ABDIc(i).ContraDist  = PathLengthToCoordinate(ABDIc(i).PreSynCoordsTransformed(ABDIc(i).isContra,:),ABDIc(i).Tree{1});
    ABDIc(i).IntegratorDist  = PathLengthToCoordinate(ABDIc(i).PreSynCoordsTransformed(ABDIc(i).isIntegrator,:),ABDIc(i).Tree{1});
    ABDIc(i).EverytingElseDist = PathLengthToCoordinate(ABDIc(i).PreSynCoordsTransformed(ABDIc(i).isEverythingElse,:),ABDIc(i).Tree{1});
end

%% plot distributions of everything

subplot(4,5,1)
h1 = histogram(vertcat(ABDr.SaccadicDist),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(1,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
subplot(4,5,6)
h1 = histogram(vertcat(ABDc.SaccadicDist),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(1,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
subplot(4,5,11)
h1 = histogram(vertcat(ABDIr.SaccadicDist),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(1,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
subplot(4,5,16)
h1 = histogram(vertcat(ABDIc.SaccadicDist),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(1,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;



subplot(4,5,2)
h1=histogram(vertcat(ABDr.VestibularDist),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(2,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
subplot(4,5,7)
h1=histogram(vertcat(ABDc.VestibularDist),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(2,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
subplot(4,5,12)
h1=histogram(vertcat(ABDIr.VestibularDist),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(2,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
subplot(4,5,17)
h1=histogram(vertcat(ABDIc.VestibularDist),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(2,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;


subplot(4,5,3)
h1=histogram(vertcat(ABDr.IntegratorDist),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(4,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
subplot(4,5,8)
h1=histogram(vertcat(ABDc.IntegratorDist),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(4,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
subplot(4,5,13)
h1=histogram(vertcat(ABDIr.IntegratorDist),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(4,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
subplot(4,5,18)
h1=histogram(vertcat(ABDIc.IntegratorDist),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(4,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;


subplot(4,5,4)
h1=histogram(vertcat(ABDr.ContraDist),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(3,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
subplot(4,5,9)
h1=histogram(vertcat(ABDc.ContraDist),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(3,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
subplot(4,5,14)
h1=histogram(vertcat(ABDIr.ContraDist),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(3,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
subplot(4,5,19)
h1=histogram(vertcat(ABDIc.ContraDist),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(3,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;


subplot(4,5,5)
h1=histogram(vertcat(ABDr.EverytingElseDist),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(5,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
subplot(4,5,10)
h1=histogram(vertcat(ABDc.EverytingElseDist),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(5,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
subplot(4,5,15)
h1=histogram(vertcat(ABDIr.EverytingElseDist),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(5,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
subplot(4,5,20)
h1=histogram(vertcat(ABDIc.EverytingElseDist),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(5,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;


%% distributions of each class onto lumped Motor classes.

subplot(4,4,1)
h1=histogram([vertcat(ABDr.SaccadicDist);vertcat(ABDc.SaccadicDist)],20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = 'm';
h1.LineWidth = 2;
set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;

h1=histogram([vertcat(ABDIr.SaccadicDist);vertcat(ABDIc.SaccadicDist)],20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = 'g';
h1.LineWidth = 2;
set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;
axis square
xlabel('Pathlength (um)');
title('Saccadic axons');

subplot(4,4,2)
h1=histogram([vertcat(ABDr.VestibularDist);vertcat(ABDc.VestibularDist)],20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = 'm';
h1.LineWidth = 2;
set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;

h1=histogram([vertcat(ABDIr.VestibularDist);vertcat(ABDIc.VestibularDist)],20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = 'g';
h1.LineWidth = 2;
set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;
axis square
title('Vestibular axons');


subplot(4,4,3)
h1=histogram([vertcat(ABDr.IntegratorDist);vertcat(ABDc.IntegratorDist)],20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = 'm';
h1.LineWidth = 2;
set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;

h1=histogram([vertcat(ABDIr.IntegratorDist);vertcat(ABDIc.IntegratorDist)],20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = 'g';
h1.LineWidth = 2;
set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;
axis square
title('Integrator axons');



subplot(4,4,4)
h1=histogram([vertcat(ABDr.ContraDist);vertcat(ABDc.ContraDist)],20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = 'm';
h1.LineWidth = 2;
set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;

h1=histogram([vertcat(ABDIr.ContraDist);vertcat(ABDIc.ContraDist)],20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = 'g';
h1.LineWidth = 2;
set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;
axis square
title('Contra axons');



subplot(4,4,5)
h1=histogram([vertcat(ABDr.EverytingElseDist);vertcat(ABDc.EverytingElseDist)],20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = 'm';
h1.LineWidth = 2;
set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;

h1=histogram([vertcat(ABDIr.EverytingElseDist);vertcat(ABDIc.EverytingElseDist)],20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = 'g';
h1.LineWidth = 2;
set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;
axis square
title('Remaining axons');

legend({'ABD','ABDi'},'Location','bestoutside');

%% Split the contributions from each of the axons classes

load('leadDiff.mat','leadDiff')
load('leadMotorDiff.mat','leadMotorDiff');
load('leadDiffAxons.mat','leadDiffAxons');

load('lagDiff.mat','lagDiff')
load('lagMotorDiff.mat','lagMotorDiff');
load('lagDiffAxons.mat','lagDiffAxons');

%%

figure;

% saccadic
subplot(4,4,1)
h1=histogram([vertcat(ABDr.SaccadicDist);vertcat(ABDc.SaccadicDist)],20,'FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = 'm';
h1.LineWidth = 2;
%set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;
h1=histogram([vertcat(ABDIr.SaccadicDist);vertcat(ABDIc.SaccadicDist)],20,'FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = 'g';
h1.LineWidth = 2;
%set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;
axis square
xlabel('Pathlength (um)');
title('Saccadic axons');

temp = [];
for i = 1:numel(ABDr_CellIDs)
    [~,ia,~] = intersect(ABDr(i).Inputs, [leadDiffAxons.Saccadic]);
    for j = 1:length(ia)
    tempdist = PathLengthToCoordinate(ABDr(i).PreSynCoordsTransformed(ia(j),:),ABDr(i).Tree{1});
    temp = [temp; tempdist];
    end
    clear ia;
    ABDr(i).LeadSaccadicDist = temp;
    temp = [];
    
    [~,ia,~] = intersect(ABDr(i).Inputs, [lagDiffAxons.Saccadic]);
    for j = 1:length(ia)
    tempdist = PathLengthToCoordinate(ABDr(i).PreSynCoordsTransformed(ia(j),:),ABDr(i).Tree{1});
    temp = [temp; tempdist];
    end
    clear ia;
    ABDr(i).LagSaccadicDist = temp;
    temp =[];
end

temp = [];
for i = 1:numel(ABDc_CellIDs)
    [~,ia,~] = intersect(ABDc(i).Inputs, [leadDiffAxons.Saccadic]);
    for j = 1:length(ia)
    tempdist = PathLengthToCoordinate(ABDc(i).PreSynCoordsTransformed(ia(j),:),ABDc(i).Tree{1});
    temp = [temp; tempdist];
    end
    clear ia;
    ABDc(i).LeadSaccadicDist = temp;
    temp = [];
    
    [~,ia,~] = intersect(ABDc(i).Inputs, [lagDiffAxons.Saccadic]);
    for j = 1:length(ia)
    tempdist = PathLengthToCoordinate(ABDc(i).PreSynCoordsTransformed(ia(j),:),ABDc(i).Tree{1});
    temp = [temp; tempdist];
    end
    clear ia;
    ABDc(i).LagSaccadicDist = temp;
    temp =[];
end
    

subplot(4,4,2)
h1=histogram([vertcat(ABDr.LeadSaccadicDist);vertcat(ABDc.LeadSaccadicDist)],20,'FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = 'r';
h1.LineWidth = 2;
%set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;
h1=histogram([vertcat(ABDr.LagSaccadicDist);vertcat(ABDc.LagSaccadicDist)],20,'FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = 'b';
h1.LineWidth = 2;
%set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;
axis square
xlabel('Pathlength (um)');
title('ABD Saccadic axons');

% ABDi Calculation

temp = [];
for i = 1:numel(ABDIr_CellIDs)
    [~,ia,~] = intersect(ABDIr(i).Inputs, [leadDiffAxons.Saccadic]);
    for j = 1:length(ia)
    tempdist = PathLengthToCoordinate(ABDIr(i).PreSynCoordsTransformed(ia(j),:),ABDIr(i).Tree{1});
    temp = [temp; tempdist];
    end
    clear ia;
    ABDIr(i).LeadSaccadicDist = temp;
    temp = [];
    
    [~,ia,~] = intersect(ABDIr(i).Inputs, [lagDiffAxons.Saccadic]);
    for j = 1:length(ia)
    tempdist = PathLengthToCoordinate(ABDIr(i).PreSynCoordsTransformed(ia(j),:),ABDIr(i).Tree{1});
    temp = [temp; tempdist];
    end
    clear ia;
    ABDIr(i).LagSaccadicDist = temp;
    temp =[];
end

temp = [];
for i = 1:numel(ABDIc_CellIDs)
    [~,ia,~] = intersect(ABDIc(i).Inputs, [leadDiffAxons.Saccadic]);
    for j = 1:length(ia)
    tempdist = PathLengthToCoordinate(ABDIc(i).PreSynCoordsTransformed(ia(j),:),ABDIc(i).Tree{1});
    temp = [temp; tempdist];
    end
    clear ia;
    ABDIc(i).LeadSaccadicDist = temp;
    temp = [];
    
    [~,ia,~] = intersect(ABDIc(i).Inputs, [lagDiffAxons.Saccadic]);
    for j = 1:length(ia)
    tempdist = PathLengthToCoordinate(ABDIc(i).PreSynCoordsTransformed(ia(j),:),ABDIc(i).Tree{1});
    temp = [temp; tempdist];
    end
    clear ia;
    ABDIc(i).LagSaccadicDist = temp;
    temp =[];
end
    

subplot(4,4,3)
h1=histogram([vertcat(ABDIr.LeadSaccadicDist);vertcat(ABDIc.LeadSaccadicDist)],20,'FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = 'r';
h1.LineWidth = 2;
%set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;
h1=histogram([vertcat(ABDIr.LagSaccadicDist);vertcat(ABDIc.LagSaccadicDist)],20,'FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = 'b';
h1.LineWidth = 2;
%set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;
axis square
xlabel('Pathlength (um)');
title('ABDi Saccadic axons');

%% Vestibular inputs distribution
subplot(4,4,5)
h1=histogram([vertcat(ABDr.VestibularDist);vertcat(ABDc.VestibularDist)],20,'FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = 'm';
h1.LineWidth = 2;
%set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;
h1=histogram([vertcat(ABDIr.VestibularDist);vertcat(ABDIc.VestibularDist)],20,'FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = 'g';
h1.LineWidth = 2;
%set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;
axis square
xlabel('Pathlength (um)');
title('Saccadic axons');

temp = [];
for i = 1:numel(ABDr_CellIDs)
    [~,ia,~] = intersect(ABDr(i).Inputs, [leadDiffAxons.Vestibular]);
    for j = 1:length(ia)
    tempdist = PathLengthToCoordinate(ABDr(i).PreSynCoordsTransformed(ia(j),:),ABDr(i).Tree{1});
    temp = [temp; tempdist];
    end
    clear ia;
    ABDr(i).LeadVestibularDist = temp;
    temp = [];
    
    [~,ia,~] = intersect(ABDr(i).Inputs, [lagDiffAxons.Vestibular]);
    for j = 1:length(ia)
    tempdist = PathLengthToCoordinate(ABDr(i).PreSynCoordsTransformed(ia(j),:),ABDr(i).Tree{1});
    temp = [temp; tempdist];
    end
    clear ia;
    ABDr(i).LagVestibularDist = temp;
    temp =[];
end

temp = [];
for i = 1:numel(ABDc_CellIDs)
    [~,ia,~] = intersect(ABDc(i).Inputs, [leadDiffAxons.Vestibular]);
    for j = 1:length(ia)
    tempdist = PathLengthToCoordinate(ABDc(i).PreSynCoordsTransformed(ia(j),:),ABDc(i).Tree{1});
    temp = [temp; tempdist];
    end
    clear ia;
    ABDc(i).LeadVestibularDist = temp;
    temp = [];
    
    [~,ia,~] = intersect(ABDc(i).Inputs, [lagDiffAxons.Vestibular]);
    for j = 1:length(ia)
    tempdist = PathLengthToCoordinate(ABDc(i).PreSynCoordsTransformed(ia(j),:),ABDc(i).Tree{1});
    temp = [temp; tempdist];
    end
    clear ia;
    ABDc(i).LagVestibularDist = temp;
    temp =[];
end
    

subplot(4,4,6)
h1=histogram([vertcat(ABDr.LeadVestibularDist);vertcat(ABDc.LeadVestibularDist)],20,'FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = 'r';
h1.LineWidth = 2;
%set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;
h1=histogram([vertcat(ABDr.LagVestibularDist);vertcat(ABDc.LagVestibularDist)],20,'FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = 'b';
h1.LineWidth = 2;
%set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;
axis square
xlabel('Pathlength (um)');
title('ABD Vestibular axons');

% ABDi Calculation

temp = [];
for i = 1:numel(ABDIr_CellIDs)
    [~,ia,~] = intersect(ABDIr(i).Inputs, [leadDiffAxons.Vestibular]);
    for j = 1:length(ia)
    tempdist = PathLengthToCoordinate(ABDIr(i).PreSynCoordsTransformed(ia(j),:),ABDIr(i).Tree{1});
    temp = [temp; tempdist];
    end
    clear ia;
    ABDIr(i).LeadVestibularDist = temp;
    temp = [];
    
    [~,ia,~] = intersect(ABDIr(i).Inputs, [lagDiffAxons.Vestibular]);
    for j = 1:length(ia)
    tempdist = PathLengthToCoordinate(ABDIr(i).PreSynCoordsTransformed(ia(j),:),ABDIr(i).Tree{1});
    temp = [temp; tempdist];
    end
    clear ia;
    ABDIr(i).LagVestibularDist = temp;
    temp =[];
end

temp = [];
for i = 1:numel(ABDIc_CellIDs)
    [~,ia,~] = intersect(ABDIc(i).Inputs, [leadDiffAxons.Vestibular]);
    for j = 1:length(ia)
    tempdist = PathLengthToCoordinate(ABDIc(i).PreSynCoordsTransformed(ia(j),:),ABDIc(i).Tree{1});
    temp = [temp; tempdist];
    end
    clear ia;
    ABDIc(i).LeadVestibularDist = temp;
    temp = [];
    
    [~,ia,~] = intersect(ABDIc(i).Inputs, [lagDiffAxons.Vestibular]);
    for j = 1:length(ia)
    tempdist = PathLengthToCoordinate(ABDIc(i).PreSynCoordsTransformed(ia(j),:),ABDIc(i).Tree{1});
    temp = [temp; tempdist];
    end
    clear ia;
    ABDIc(i).LagVestibularDist = temp;
    temp =[];
end
    

subplot(4,4,7)
h1=histogram([vertcat(ABDIr.LeadVestibularDist);vertcat(ABDIc.LeadVestibularDist)],20,'FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = 'r';
h1.LineWidth = 2;
%set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;
h1=histogram([vertcat(ABDIr.LagVestibularDist);vertcat(ABDIc.LagVestibularDist)],20,'FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = 'b';
h1.LineWidth = 2;
%set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;
axis square
xlabel('Pathlength (um)');
title('ABDi Vestibular axons');

%% Integrator axons distribution
subplot(4,4,9)
h1=histogram([vertcat(ABDr.IntegratorDist);vertcat(ABDc.IntegratorDist)],20,'FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = 'm';
h1.LineWidth = 2;
%set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;
h1=histogram([vertcat(ABDIr.IntegratorDist);vertcat(ABDIc.IntegratorDist)],20,'FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = 'g';
h1.LineWidth = 2;
%set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;
axis square
xlabel('Pathlength (um)');
title('Integrator axons');

temp = [];
for i = 1:numel(ABDr_CellIDs)
    [~,ia,~] = intersect(ABDr(i).Inputs, [leadDiffAxons.Integrator]);
    for j = 1:length(ia)
    tempdist = PathLengthToCoordinate(ABDr(i).PreSynCoordsTransformed(ia(j),:),ABDr(i).Tree{1});
    temp = [temp; tempdist];
    end
    clear ia;
    ABDr(i).LeadIntegratorDist = temp;
    temp = [];
    
    [~,ia,~] = intersect(ABDr(i).Inputs, [lagDiffAxons.Integrator]);
    for j = 1:length(ia)
    tempdist = PathLengthToCoordinate(ABDr(i).PreSynCoordsTransformed(ia(j),:),ABDr(i).Tree{1});
    temp = [temp; tempdist];
    end
    clear ia;
    ABDr(i).LagIntegratorDist = temp;
    temp =[];
end

temp = [];
for i = 1:numel(ABDc_CellIDs)
    [~,ia,~] = intersect(ABDc(i).Inputs, [leadDiffAxons.Integrator]);
    for j = 1:length(ia)
    tempdist = PathLengthToCoordinate(ABDc(i).PreSynCoordsTransformed(ia(j),:),ABDc(i).Tree{1});
    temp = [temp; tempdist];
    end
    clear ia;
    ABDc(i).LeadIntegratorDist = temp;
    temp = [];
    
    [~,ia,~] = intersect(ABDc(i).Inputs, [lagDiffAxons.Integrator]);
    for j = 1:length(ia)
    tempdist = PathLengthToCoordinate(ABDc(i).PreSynCoordsTransformed(ia(j),:),ABDc(i).Tree{1});
    temp = [temp; tempdist];
    end
    clear ia;
    ABDc(i).LagIntegratorDist = temp;
    temp =[];
end
    

subplot(4,4,10)
h1=histogram([vertcat(ABDr.LeadIntegratorDist);vertcat(ABDc.LeadIntegratorDist)],20,'FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = 'r';
h1.LineWidth = 2;
%set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;
h1=histogram([vertcat(ABDr.LagIntegratorDist);vertcat(ABDc.LagIntegratorDist)],20,'FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = 'b';
h1.LineWidth = 2;
%set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;
axis square
xlabel('Pathlength (um)');
title('ABD Integrator axons');

% ABDi Calculation

temp = [];
for i = 1:numel(ABDIr_CellIDs)
    [~,ia,~] = intersect(ABDIr(i).Inputs, [leadDiffAxons.Integrator]);
    for j = 1:length(ia)
    tempdist = PathLengthToCoordinate(ABDIr(i).PreSynCoordsTransformed(ia(j),:),ABDIr(i).Tree{1});
    temp = [temp; tempdist];
    end
    clear ia;
    ABDIr(i).LeadIntegratorDist = temp;
    temp = [];
    
    [~,ia,~] = intersect(ABDIr(i).Inputs, [lagDiffAxons.Integrator]);
    for j = 1:length(ia)
    tempdist = PathLengthToCoordinate(ABDIr(i).PreSynCoordsTransformed(ia(j),:),ABDIr(i).Tree{1});
    temp = [temp; tempdist];
    end
    clear ia;
    ABDIr(i).LagIntegratorDist = temp;
    temp =[];
end

temp = [];
for i = 1:numel(ABDIc_CellIDs)
    [~,ia,~] = intersect(ABDIc(i).Inputs, [leadDiffAxons.Integrator]);
    for j = 1:length(ia)
    tempdist = PathLengthToCoordinate(ABDIc(i).PreSynCoordsTransformed(ia(j),:),ABDIc(i).Tree{1});
    temp = [temp; tempdist];
    end
    clear ia;
    ABDIc(i).LeadIntegratorDist = temp;
    temp = [];
    
    [~,ia,~] = intersect(ABDIc(i).Inputs, [lagDiffAxons.Integrator]);
    for j = 1:length(ia)
    tempdist = PathLengthToCoordinate(ABDIc(i).PreSynCoordsTransformed(ia(j),:),ABDIc(i).Tree{1});
    temp = [temp; tempdist];
    end
    clear ia;
    ABDIc(i).LagIntegratorDist = temp;
    temp =[];
end
    

subplot(4,4,11)
h1=histogram([vertcat(ABDIr.LeadIntegratorDist);vertcat(ABDIc.LeadIntegratorDist)],20,'FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = 'r';
h1.LineWidth = 2;
%set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;
h1=histogram([vertcat(ABDIr.LagIntegratorDist);vertcat(ABDIc.LagIntegratorDist)],20,'FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = 'b';
h1.LineWidth = 2;
%set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;
axis square
xlabel('Pathlength (um)');
title('ABDi Integrator axons');

%% Contra axons distribution
subplot(4,4,13)
h1=histogram([vertcat(ABDr.ContraDist);vertcat(ABDc.ContraDist)],20,'FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = 'm';
h1.LineWidth = 2;
%set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;
h1=histogram([vertcat(ABDIr.ContraDist);vertcat(ABDIc.ContraDist)],20,'FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = 'g';
h1.LineWidth = 2;
%set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;
axis square
xlabel('Pathlength (um)');
title('Contra axons');

temp = [];
for i = 1:numel(ABDr_CellIDs)
    [~,ia,~] = intersect(ABDr(i).Inputs, [leadDiffAxons.Contra]);
    for j = 1:length(ia)
    tempdist = PathLengthToCoordinate(ABDr(i).PreSynCoordsTransformed(ia(j),:),ABDr(i).Tree{1});
    temp = [temp; tempdist];
    end
    clear ia;
    ABDr(i).LeadContraDist = temp;
    temp = [];
    
    [~,ia,~] = intersect(ABDr(i).Inputs, [lagDiffAxons.Contra]);
    for j = 1:length(ia)
    tempdist = PathLengthToCoordinate(ABDr(i).PreSynCoordsTransformed(ia(j),:),ABDr(i).Tree{1});
    temp = [temp; tempdist];
    end
    clear ia;
    ABDr(i).LagContraDist = temp;
    temp =[];
end

temp = [];
for i = 1:numel(ABDc_CellIDs)
    [~,ia,~] = intersect(ABDc(i).Inputs, [leadDiffAxons.Contra]);
    for j = 1:length(ia)
    tempdist = PathLengthToCoordinate(ABDc(i).PreSynCoordsTransformed(ia(j),:),ABDc(i).Tree{1});
    temp = [temp; tempdist];
    end
    clear ia;
    ABDc(i).LeadContraDist = temp;
    temp = [];
    
    [~,ia,~] = intersect(ABDc(i).Inputs, [lagDiffAxons.Contra]);
    for j = 1:length(ia)
    tempdist = PathLengthToCoordinate(ABDc(i).PreSynCoordsTransformed(ia(j),:),ABDc(i).Tree{1});
    temp = [temp; tempdist];
    end
    clear ia;
    ABDc(i).LagContraDist = temp;
    temp =[];
end
    

subplot(4,4,14)
h1=histogram([vertcat(ABDr.LeadContraDist);vertcat(ABDc.LeadContraDist)],20,'FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = 'r';
h1.LineWidth = 2;
%set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;
h1=histogram([vertcat(ABDr.LagContraDist);vertcat(ABDc.LagContraDist)],20,'FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = 'b';
h1.LineWidth = 2;
%set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;
axis square
xlabel('Pathlength (um)');
title('ABD Contra axons');

% ABDi Calculation

temp = [];
for i = 1:numel(ABDIr_CellIDs)
    [~,ia,~] = intersect(ABDIr(i).Inputs, [leadDiffAxons.Contra]);
    for j = 1:length(ia)
    tempdist = PathLengthToCoordinate(ABDIr(i).PreSynCoordsTransformed(ia(j),:),ABDIr(i).Tree{1});
    temp = [temp; tempdist];
    end
    clear ia;
    ABDIr(i).LeadContraDist = temp;
    temp = [];
    
    [~,ia,~] = intersect(ABDIr(i).Inputs, [lagDiffAxons.Contra]);
    for j = 1:length(ia)
    tempdist = PathLengthToCoordinate(ABDIr(i).PreSynCoordsTransformed(ia(j),:),ABDIr(i).Tree{1});
    temp = [temp; tempdist];
    end
    clear ia;
    ABDIr(i).LagContraDist = temp;
    temp =[];
end

temp = [];
for i = 1:numel(ABDIc_CellIDs)
    [~,ia,~] = intersect(ABDIc(i).Inputs, [leadDiffAxons.Contra]);
    for j = 1:length(ia)
    tempdist = PathLengthToCoordinate(ABDIc(i).PreSynCoordsTransformed(ia(j),:),ABDIc(i).Tree{1});
    temp = [temp; tempdist];
    end
    clear ia;
    ABDIc(i).LeadContraDist = temp;
    temp = [];
    
    [~,ia,~] = intersect(ABDIc(i).Inputs, [lagDiffAxons.Contra]);
    for j = 1:length(ia)
    tempdist = PathLengthToCoordinate(ABDIc(i).PreSynCoordsTransformed(ia(j),:),ABDIc(i).Tree{1});
    temp = [temp; tempdist];
    end
    clear ia;
    ABDIc(i).LagContraDist = temp;
    temp =[];
end
    

subplot(4,4,15)
h1=histogram([vertcat(ABDIr.LeadContraDist);vertcat(ABDIc.LeadContraDist)],20,'FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = 'r';
h1.LineWidth = 2;
%set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;
h1=histogram([vertcat(ABDIr.LagContraDist);vertcat(ABDIc.LagContraDist)],20,'FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = 'b';
h1.LineWidth = 2;
%set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;
axis square
xlabel('Pathlength (um)');
title('ABDi Contra axons');

%% Remaining axon distribution

subplot(4,4,1)
h1=histogram([vertcat(ABDr.EverythingElseDist);vertcat(ABDc.EverythingElseDist)],20,'FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = 'm';
h1.LineWidth = 2;
%set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;
h1=histogram([vertcat(ABDIr.EverythingElseDist);vertcat(ABDIc.EverythingElseDist)],20,'FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = 'g';
h1.LineWidth = 2;
%set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;
axis square
xlabel('Pathlength (um)');
title('Remaining axons');

temp = [];
for i = 1:numel(ABDr_CellIDs)
    [~,ia,~] = intersect(ABDr(i).Inputs, [leadDiffAxons.EverythingElse]);
    for j = 1:length(ia)
    tempdist = PathLengthToCoordinate(ABDr(i).PreSynCoordsTransformed(ia(j),:),ABDr(i).Tree{1});
    temp = [temp; tempdist];
    end
    clear ia;
    ABDr(i).LeadEverythingElseDist = temp;
    temp = [];
    
    [~,ia,~] = intersect(ABDr(i).Inputs, [lagDiffAxons.EverythingElse]);
    for j = 1:length(ia)
    tempdist = PathLengthToCoordinate(ABDr(i).PreSynCoordsTransformed(ia(j),:),ABDr(i).Tree{1});
    temp = [temp; tempdist];
    end
    clear ia;
    ABDr(i).LagEverythingElseDist = temp;
    temp =[];
end

temp = [];
for i = 1:numel(ABDc_CellIDs)
    [~,ia,~] = intersect(ABDc(i).Inputs, [leadDiffAxons.EverythingElse]);
    for j = 1:length(ia)
    tempdist = PathLengthToCoordinate(ABDc(i).PreSynCoordsTransformed(ia(j),:),ABDc(i).Tree{1});
    temp = [temp; tempdist];
    end
    clear ia;
    ABDc(i).LeadEverythingElseDist = temp;
    temp = [];
    
    [~,ia,~] = intersect(ABDc(i).Inputs, [lagDiffAxons.EverythingElse]);
    for j = 1:length(ia)
    tempdist = PathLengthToCoordinate(ABDc(i).PreSynCoordsTransformed(ia(j),:),ABDc(i).Tree{1});
    temp = [temp; tempdist];
    end
    clear ia;
    ABDc(i).LagEverythingElseDist = temp;
    temp =[];
end
    

subplot(4,4,2)
h1=histogram([vertcat(ABDr.LeadEverythingElseDist);vertcat(ABDc.LeadEverythingElseDist)],20,'FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = 'r';
h1.LineWidth = 2;
%set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;
h1=histogram([vertcat(ABDr.LagEverythingElseDist);vertcat(ABDc.LagEverythingElseDist)],20,'FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = 'b';
h1.LineWidth = 2;
%set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;
axis square
xlabel('Pathlength (um)');
title('ABD EverythingElse axons');

% ABDi Calculation

temp = [];
for i = 1:numel(ABDIr_CellIDs)
    [~,ia,~] = intersect(ABDIr(i).Inputs, [leadDiffAxons.EverythingElse]);
    for j = 1:length(ia)
    tempdist = PathLengthToCoordinate(ABDIr(i).PreSynCoordsTransformed(ia(j),:),ABDIr(i).Tree{1});
    temp = [temp; tempdist];
    end
    clear ia;
    ABDIr(i).LeadEverythingElseDist = temp;
    temp = [];
    
    [~,ia,~] = intersect(ABDIr(i).Inputs, [lagDiffAxons.EverythingElse]);
    for j = 1:length(ia)
    tempdist = PathLengthToCoordinate(ABDIr(i).PreSynCoordsTransformed(ia(j),:),ABDIr(i).Tree{1});
    temp = [temp; tempdist];
    end
    clear ia;
    ABDIr(i).LagEverythingElseDist = temp;
    temp =[];
end

temp = [];
for i = 1:numel(ABDIc_CellIDs)
    [~,ia,~] = intersect(ABDIc(i).Inputs, [leadDiffAxons.EverythingElse]);
    for j = 1:length(ia)
    tempdist = PathLengthToCoordinate(ABDIc(i).PreSynCoordsTransformed(ia(j),:),ABDIc(i).Tree{1});
    temp = [temp; tempdist];
    end
    clear ia;
    ABDIc(i).LeadEverythingElseDist = temp;
    temp = [];
    
    [~,ia,~] = intersect(ABDIc(i).Inputs, [lagDiffAxons.EverythingElse]);
    for j = 1:length(ia)
    tempdist = PathLengthToCoordinate(ABDIc(i).PreSynCoordsTransformed(ia(j),:),ABDIc(i).Tree{1});
    temp = [temp; tempdist];
    end
    clear ia;
    ABDIc(i).LagEverythingElseDist = temp;
    temp =[];
end
    

subplot(4,4,3)
h1=histogram([vertcat(ABDIr.LeadEverythingElseDist);vertcat(ABDIc.LeadEverythingElseDist)],20,'FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = 'r';
h1.LineWidth = 2;
%set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;
h1=histogram([vertcat(ABDIr.LagEverythingElseDist);vertcat(ABDIc.LagEverythingElseDist)],20,'FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = 'b';
h1.LineWidth = 2;
%set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;
axis square
xlabel('Pathlength (um)');
title('ABDi EverythingElse axons');

%% plot location on individual Abducens neurons

figure;
for i = 1:numel(ABDr_CellIDs)
    subplot(4,4,i)
plot_tree(ABDr(i).Tree{1},'k',[],[],[],'-3l');
hold on;
scatter3(ABDr(i).PreSynCoordsTransformed(ABDr(i).isSaccadic,1), ABDr(i).PreSynCoordsTransformed(ABDr(i).isSaccadic,2), ABDr(i).PreSynCoordsTransformed(ABDr(i).isSaccadic,3),'o',...
    'MarkerFaceColor',colors(1,:),'MarkerEdgeColor','none')
scatter3(ABDr(i).PreSynCoordsTransformed(ABDr(i).isContra,1), ABDr(i).PreSynCoordsTransformed(ABDr(i).isContra,2), ABDr(i).PreSynCoordsTransformed(ABDr(i).isContra,3),'o', ...
     'MarkerFaceColor',colors(2,:),'MarkerEdgeColor','none')
scatter3(ABDr(i).PreSynCoordsTransformed(ABDr(i).isIntegrator,1), ABDr(i).PreSynCoordsTransformed(ABDr(i).isIntegrator,2), ABDr(i).PreSynCoordsTransformed(ABDr(i).isIntegrator,3),'o', ...
 'MarkerFaceColor',colors(3,:),'MarkerEdgeColor','none')
title(ABDr_CellIDs(i));
axis off;
end
sgtitle('ABDr');

figure;
for i = 1:numel(ABDc_CellIDs)
    subplot(4,4,i)
plot_tree(ABDc(i).Tree{1},'k',[],[],[],'-3l');
hold on;
scatter3(ABDc(i).PreSynCoordsTransformed(ABDc(i).isSaccadic,1), ABDc(i).PreSynCoordsTransformed(ABDc(i).isSaccadic,2), ABDc(i).PreSynCoordsTransformed(ABDc(i).isSaccadic,3),'o',...
    'MarkerFaceColor',colors(1,:),'MarkerEdgeColor','none')
scatter3(ABDc(i).PreSynCoordsTransformed(ABDc(i).isContra,1), ABDc(i).PreSynCoordsTransformed(ABDc(i).isContra,2), ABDc(i).PreSynCoordsTransformed(ABDc(i).isContra,3),'o', ...
     'MarkerFaceColor',colors(2,:),'MarkerEdgeColor','none')
scatter3(ABDc(i).PreSynCoordsTransformed(ABDc(i).isIntegrator,1), ABDc(i).PreSynCoordsTransformed(ABDc(i).isIntegrator,2), ABDc(i).PreSynCoordsTransformed(ABDc(i).isIntegrator,3),'o', ...
 'MarkerFaceColor',colors(3,:),'MarkerEdgeColor','none')
title(ABDc_CellIDs(i));
axis off;
end
sgtitle('ABDc');

figure;
for i = 1:numel(ABDIr_CellIDs)
    subplot(4,4,i)
plot_tree(ABDIr(i).Tree{1},'k',[],[],[],'-3l');
hold on;
scatter3(ABDIr(i).PreSynCoordsTransformed(ABDIr(i).isSaccadic,1), ABDIr(i).PreSynCoordsTransformed(ABDIr(i).isSaccadic,2), ABDIr(i).PreSynCoordsTransformed(ABDIr(i).isSaccadic,3),'o',...
    'MarkerFaceColor',colors(1,:),'MarkerEdgeColor','none')
scatter3(ABDIr(i).PreSynCoordsTransformed(ABDIr(i).isContra,1), ABDIr(i).PreSynCoordsTransformed(ABDIr(i).isContra,2), ABDIr(i).PreSynCoordsTransformed(ABDIr(i).isContra,3),'o', ...
     'MarkerFaceColor',colors(2,:),'MarkerEdgeColor','none')
scatter3(ABDIr(i).PreSynCoordsTransformed(ABDIr(i).isIntegrator,1), ABDIr(i).PreSynCoordsTransformed(ABDIr(i).isIntegrator,2), ABDIr(i).PreSynCoordsTransformed(ABDIr(i).isIntegrator,3),'o', ...
 'MarkerFaceColor',colors(3,:),'MarkerEdgeColor','none')
title(ABDIr_CellIDs(i));
axis off;
end
sgtitle('ABDIr');

figure;
for i = 1:numel(ABDIc_CellIDs)
    subplot(4,4,i)
plot_tree(ABDIc(i).Tree{1},'k',[],[],[],'-3l');
hold on;
scatter3(ABDIc(i).PreSynCoordsTransformed(ABDIc(i).isSaccadic,1), ABDIc(i).PreSynCoordsTransformed(ABDIc(i).isSaccadic,2), ABDIc(i).PreSynCoordsTransformed(ABDIc(i).isSaccadic,3),'o',...
    'MarkerFaceColor',colors(1,:),'MarkerEdgeColor','none')
scatter3(ABDIc(i).PreSynCoordsTransformed(ABDIc(i).isContra,1), ABDIc(i).PreSynCoordsTransformed(ABDIc(i).isContra,2), ABDIc(i).PreSynCoordsTransformed(ABDIc(i).isContra,3),'o', ...
     'MarkerFaceColor',colors(2,:),'MarkerEdgeColor','none')
scatter3(ABDIc(i).PreSynCoordsTransformed(ABDIc(i).isIntegrator,1), ABDIc(i).PreSynCoordsTransformed(ABDIc(i).isIntegrator,2), ABDIc(i).PreSynCoordsTransformed(ABDIc(i).isIntegrator,3),'o', ...
 'MarkerFaceColor',colors(3,:),'MarkerEdgeColor','none')
title(ABDIc_CellIDs(i));
axis off;
end
sgtitle('ABDIc');


%% 
% complete recontrustucon for 
% ABDr - 77710, 77648
% ABDIr - 78556
% ABDIc - 77641, 77148

% locate the neuron of interest
l = find(ABDr_CellIDs == 77710);
i = l;

PreSyn77710Saccade =  PathLengthToCoordinate(ABDr(i).PreSynCoordsTransformed(ABDr(i).isSaccadic,:),ABDr(i).Tree{1});
PreSyn77710Vestibular =  PathLengthToCoordinate(ABDr(i).PreSynCoordsTransformed(ABDr(i).isVestibular,:),ABDr(i).Tree{1});
PreSyn77710Contra = PathLengthToCoordinate(ABDr(i).PreSynCoordsTransformed(ABDr(i).isContra,:),ABDr(i).Tree{1});
PreSyn77710Integrator =  PathLengthToCoordinate(ABDr(i).PreSynCoordsTransformed(ABDr(i).isIntegrator,:),ABDr(i).Tree{1});
PreSyn77710EverythginElse = PathLengthToCoordinate(ABDr(i).PreSynCoordsTransformed(ABDr(i).isEverythingElse,:),ABDr(i).Tree{1});

figure()

h1 = histogram(PreSyn77710Saccade,20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(1,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;

hold on;
h1 = histogram(PreSyn77710Vestibular,20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(2,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;

h1 = histogram(PreSyn77710Contra,20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(3,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;

h1 = histogram(PreSyn77710Integrator,20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(4,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;

h1 = histogram(PreSyn77710EverythginElse,20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(5,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;







