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

for i = 1:numel(ABDr_CellIDs)
    ABDr(i) = InputsByClass(ABDr_CellIDs(i),df);  
end

for i = 1:numel(ABDc_CellIDs)
    ABDc(i) = InputsByClass(ABDc_CellIDs(i),df);
end

for i = 1:numel(ABDIr_CellIDs)
    ABDIr(i) = InputsByClass(ABDIr_CellIDs(i),df);
end

for i = 1:numel(ABDIc_CellIDs)
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
figure();
ABDrSaccadicCoords = [];
ABDrVestibularCoords = [];
ABDrContraCoords = [];
ABDrIntegratorCoords = [];
ABDrEverythingElseCoords = [];

for i = 1:numel(ABDr_CellIDs)
    ABDrSaccadicCoords  = [ABDrSaccadicCoords; TransformPoints(ABDr(i).PreSynCoords(ABDr(i).isSaccadic,:),0)-ABDr(i).Origin ];
    ABDrVestibularCoords  = [ABDrVestibularCoords; TransformPoints(ABDr(i).PreSynCoords(ABDr(i).isVestibular,:),0)-ABDr(i).Origin];
    ABDrContraCoords  = [ABDrContraCoords; TransformPoints(ABDr(i).PreSynCoords(ABDr(i).isContra,:),0)-ABDr(i).Origin];
    ABDrIntegratorCoords  = [ABDrIntegratorCoords; TransformPoints(ABDr(i).PreSynCoords(ABDr(i).isIntegrator,:),0)-ABDr(i).Origin]; 
    ABDrEverythingElseCoords = [ABDrEverythingElseCoords;TransformPoints(ABDr(i).PreSynCoords(ABDr(i).isEverythingElse,:),0)-ABDr(i).Origin];
end

ABDrSaccadicEucDist = (pdist2(zeros(size(ABDrSaccadicCoords,1),3),ABDrSaccadicCoords));
ABDrVestibularEucDist = (pdist2(zeros(size(ABDrVestibularCoords,1),3),ABDrVestibularCoords));
ABDrContraEucDist = (pdist2(zeros(size(ABDrContraCoords,1),3),ABDrContraCoords));
ABDrIntegratorEucDist = (pdist2(zeros(size(ABDrIntegratorCoords,1),3),ABDrIntegratorCoords));
ABDrEverytingElseEucDist = (pdist2(zeros(size(ABDrEverythingElseCoords,1),3),ABDrEverythingElseCoords));

% ABDc

ABDcSaccadicCoords = [];
ABDcVestibularCoords = [];
ABDcContraCoords = [];
ABDcIntegratorCoords = [];
ABDcEverythingElseCoords = [];

for i = 1:numel(ABDc_CellIDs)
    ABDcSaccadicCoords  = [ABDcSaccadicCoords; TransformPoints(ABDc(i).PreSynCoords(ABDc(i).isSaccadic,:),0)-ABDc(i).Origin];
    ABDcVestibularCoords  = [ABDcVestibularCoords; TransformPoints(ABDc(i).PreSynCoords(ABDc(i).isVestibular,:),0)-ABDc(i).Origin];
    ABDcContraCoords  = [ABDcContraCoords;TransformPoints( ABDc(i).PreSynCoords(ABDc(i).isContra,:),0)-ABDc(i).Origin];
    ABDcIntegratorCoords  = [ABDcIntegratorCoords; TransformPoints(ABDc(i).PreSynCoords(ABDc(i).isIntegrator,:),0)-ABDc(i).Origin];
    ABDcEverythingElseCoords = [ABDcEverythingElseCoords;TransformPoints(ABDc(i).PreSynCoords(ABDc(i).isEverythingElse,:),0)-ABDc(i).Origin];
end

ABDcSaccadicEucDist = (pdist2(zeros(size(ABDcSaccadicCoords,1),3),ABDcSaccadicCoords));
ABDcVestibularEucDist = (pdist2(zeros(size(ABDcVestibularCoords,1),3),ABDcVestibularCoords));
ABDcContraEucDist = (pdist2(zeros(size(ABDcContraCoords,1),3),ABDcContraCoords));
ABDcIntegratorEucDist = (pdist2(zeros(size(ABDcIntegratorCoords,1),3),ABDcIntegratorCoords));
ABDcEverytingElseEucDist = (pdist2(zeros(size(ABDcEverythingElseCoords,1),3),ABDcEverythingElseCoords));


% ABDIr

ABDIrSaccadicCoords = [];
ABDIrVestibularCoords = [];
ABDIrContraCoords = [];
ABDIrIntegratorCoords = [];
ABDIrEverythingElseCoords = [];

for i = 1:numel(ABDIr_CellIDs)
    ABDIrSaccadicCoords  = [ABDIrSaccadicCoords; TransformPoints(ABDIr(i).PreSynCoords(ABDIr(i).isSaccadic,:),0)-ABDIr(i).Origin];
    ABDIrVestibularCoords  = [ABDIrVestibularCoords; TransformPoints(ABDIr(i).PreSynCoords(ABDIr(i).isVestibular,:),0)-ABDIr(i).Origin];
    ABDIrContraCoords  = [ABDIrContraCoords; TransformPoints(ABDIr(i).PreSynCoords(ABDIr(i).isContra,:),0)-ABDIr(i).Origin];
    ABDIrIntegratorCoords  = [ABDIrIntegratorCoords; TransformPoints(ABDIr(i).PreSynCoords(ABDIr(i).isIntegrator,:),0)-ABDIr(i).Origin];  
    ABDIrEverythingElseCoords  = [ABDIrEverythingElseCoords;TransformPoints(ABDIr(i).PreSynCoords(ABDIr(i).isEverythingElse,:),0)-ABDIr(i).Origin];      
end

ABDIrSaccadicEucDist = (pdist2(zeros(size(ABDIrSaccadicCoords,1),3),ABDIrSaccadicCoords));
ABDIrVestibularEucDist = (pdist2(zeros(size(ABDIrVestibularCoords,1),3),ABDIrVestibularCoords));
ABDIrContraEucDist = (pdist2(zeros(size(ABDIrContraCoords,1),3),ABDIrContraCoords));
ABDIrIntegratorEucDist = (pdist2(zeros(size(ABDIrIntegratorCoords,1),3),ABDIrIntegratorCoords));
ABDIrEverytingElseEucDist = (pdist2(zeros(size(ABDIrEverythingElseCoords,1),3),ABDIrEverythingElseCoords));


% ABDIc


ABDIcSaccadicCoords = [];
ABDIcVestibularCoords = [];
ABDIcContraCoords = [];
ABDIcIntegratorCoords = [];
ABDIcEverythingElseCoords = [];

for i = 1:numel(ABDIc_CellIDs)
    ABDIcSaccadicCoords  = [ABDIcSaccadicCoords; TransformPoints(ABDIc(i).PreSynCoords(ABDIc(i).isSaccadic,:),0) - ABDIc(i).Origin];
    ABDIcVestibularCoords  = [ABDIcVestibularCoords; TransformPoints(ABDIc(i).PreSynCoords(ABDIc(i).isVestibular,:),0) - ABDIc(i).Origin];
    ABDIcContraCoords  = [ABDIcContraCoords; TransformPoints(ABDIc(i).PreSynCoords(ABDIc(i).isContra,:),0) - ABDIc(i).Origin];
    ABDIcIntegratorCoords  = [ABDIcIntegratorCoords; TransformPoints(ABDIc(i).PreSynCoords(ABDIc(i).isIntegrator,:),0) - ABDIc(i).Origin];  
    ABDIcEverythingElseCoords  = [ABDIcEverythingElseCoords; TransformPoints(ABDIc(i).PreSynCoords(ABDIc(i).isEverythingElse,:),0) - ABDIc(i).Origin];  

end

ABDIcSaccadicEucDist = (pdist2(zeros(size(ABDIcSaccadicCoords,1),3),ABDIcSaccadicCoords));
ABDIcVestibularEucDist = (pdist2(zeros(size(ABDIcVestibularCoords,1),3),ABDIcVestibularCoords));
ABDIcContraEucDist = (pdist2(zeros(size(ABDIcContraCoords,1),3),ABDIcContraCoords));
ABDIcIntegratorEucDist = (pdist2(zeros(size(ABDIcIntegratorCoords,1),3),ABDIcIntegratorCoords));
ABDIcEverytingElseEucDist = (pdist2(zeros(size(ABDIcEverythingElseCoords,1),3),ABDIcEverythingElseCoords));

% plot distributions of everything

subplot(4,5,1)
h1 = histogram(ABDrSaccadicEucDist(1,:),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(1,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,60]);
box off;
subplot(4,5,6)
h1 = histogram(ABDcSaccadicEucDist(1,:),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(1,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,60]);
box off;
subplot(4,5,11)
h1 = histogram(ABDIrSaccadicEucDist(1,:),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(1,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,60]);
box off;
subplot(4,5,16)
h1 = histogram(ABDIcSaccadicEucDist(1,:),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(1,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,60]);
box off;


subplot(4,5,2)
h1=histogram(ABDrVestibularEucDist(1,:),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(2,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,60]);
box off;
subplot(4,5,7)
h1=histogram(ABDcVestibularEucDist(1,:),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(2,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,60]);
box off;
subplot(4,5,12)
h1=histogram(ABDIrVestibularEucDist(1,:),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(2,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,60]);
box off;
subplot(4,5,17)
h1=histogram(ABDIcVestibularEucDist(1,:),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(2,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,60]);
box off;


subplot(4,5,3)
h1=histogram(ABDrIntegratorEucDist(1,:),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(4,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,60]);
box off;
subplot(4,5,8)
h1=histogram(ABDcIntegratorEucDist(1,:),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(4,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,60]);
box off;
subplot(4,5,13)
h1=histogram(ABDIrIntegratorEucDist(1,:),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(4,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,60]);
box off;
subplot(4,5,18)
h1=histogram(ABDIcIntegratorEucDist(1,:),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(4,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,60]);
box off;


subplot(4,5,4)
h1=histogram(ABDrContraEucDist(1,:),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(3,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,60]);
box off;
subplot(4,5,9)
h1=histogram(ABDcContraEucDist(1,:),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(3,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,60]);
box off;
subplot(4,5,14)
h1=histogram(ABDIrContraEucDist(1,:),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(3,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,60]);
box off;
subplot(4,5,19)
h1=histogram(ABDIcContraEucDist(1,:),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(3,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,60]);
box off;


subplot(4,5,5)
h1=histogram(ABDrEverytingElseEucDist(1,:),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(5,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,60]);
box off;
subplot(4,5,10)
h1=histogram(ABDcEverytingElseEucDist(1,:),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(5,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,60]);
box off;
subplot(4,5,15)
h1=histogram(ABDIrEverytingElseEucDist(1,:),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(5,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,60]);
box off;
subplot(4,5,20)
h1=histogram(ABDIcEverytingElseEucDist(1,:),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(5,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,60]);
box off;


%% 
% complete recontrustucon for 
% ABDr - 77710, 77648
% ABDIr - 78556
% ABDIc - 77641, 77148

% locate the neuron of interest
l = find(ABDr_CellIDs == 77710);
i = l;

PreSyn77710Saccade =  TransformPoints(ABDIc(i).PreSynCoords(ABDIc(i).isSaccadic,:),0)-ABDIc(i).Origin;
PreSyn77710Vestibular =  TransformPoints(ABDIc(i).PreSynCoords(ABDIc(i).isVestibular,:),0)-ABDIc(i).Origin;
PreSyn77710Contra =  TransformPoints(ABDIc(i).PreSynCoords(ABDIc(i).isContra,:),0)-ABDIc(i).Origin;
PreSyn77710Integrator =  TransformPoints(ABDIc(i).PreSynCoords(ABDIc(i).isIntegrator,:),0)-ABDIc(i).Origin;
PreSyn77710EverythginElse = TransformPoints(ABDIc(i).PreSynCoords(ABDIc(i).isEverythingElse,:),0)-ABDIc(i).Origin;


figure()

h1=histfit(PreSyn77710Saccade(:,1),50,'kernel')
h1(1).FaceColor = 'none';
h1(1).EdgeColor = 'none';
h1(2).Color = colors(1,:);
hold on;
h1=histfit(PreSyn77710Vestibular(:,1),50,'kernel');
h1(1).FaceColor = 'none';
h1(1).EdgeColor = 'none';
h1(2).Color = colors(2,:);
h1=histfit(PreSyn77710Contra(:,1),50,'kernel');
h1(1).FaceColor = 'none';
h1(1).EdgeColor = 'none';
h1(2).Color = colors(3,:);
h1=histfit(PreSyn77710Integrator(:,1),50,'kernel');
h1(1).FaceColor = 'none';
h1(1).EdgeColor = 'none';
h1(2).Color = colors(4,:);
h1=histfit(PreSyn77710EverythginElse(:,1),50,'kernel');
h1(1).FaceColor = 'none';
h1(1).EdgeColor = 'none';
h1(2).Color = colors(5,:);

box off;




