clc;
clear;

if ismac
    addpath(genpath('/Users/ashwin/Documents/'));
    df = readtable('//Users/ashwin/Google Drive/Zfish/SynapseDetector/04152019.csv');
    fname  = '/Users/ashwin/Google Drive/Zfish/LowEMtoHighEM/SWC_all/consensus-20180920/swc/';

else
    addpath(genpath('/usr/people/ashwinv/seungmount/research/Ashwin/Scripts'));
    df = readtable('/usr/people/ashwinv/seungmount/research/Ashwin/SynapseDetector/09202018.csv');
    %ml = readtable('/usr/people/ashwinv/seungmount/research/Ashwin/SynapseDetector/ManualSortList-102218.csv');
    fname =  '/usr/people/ashwinv/seungmount/research/Ashwin/Z-Brain/LowEMtoHighEM/SWC_all/consensus-20180920/swc/'

end

colors = cbrewer('qual','Dark2',10);   

temp1 = colorcet('CBTL1','N',5); %  linear-tritanopic_krjcw_5-98_c46_n256
temp2 = colorcet('CBL2','N',5);  %  linear-protanopic-deuteranopic_kbw_5-98_c40_n256

leadColor = temp1(3,:); % dark red, CB safe
lagColor = temp2(3,:); % dark blue, CB safe

ABDr_CellIDs = [82140,82145,82143,77648,82146,77302,77301, 82194,77705,77710,77305,77672,77300,77709];
ABDc_CellIDs = [82213,77654,77646,82212,77295,81172,77657,82197,77682,77154, 77652,77658,77628,77292,77688,82195,77296];
ABDIr_CellIDs = [78553,77150,77886,78547,77631,77158,77665,78552,77668,77618,77634];
ABDIc_CellIDs = [78574,79051,77148,77625,77144,77692,77641,79066,77640,77643];

% complete recontrustucon for 
% ABDr - 77710, 77648
% ABDIr - 78556
% ABDIc - 77641, 77148
%
load ABDPutativeSaccadic.mat
load ABDiPutativeSaccadic.mat

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

    
    save('ABDr.mat','ABDr');
    save('ABDc.mat','ABDc');
    save('ABDIr.mat','ABDIr');
    save('ABDIc.mat','ABDIc');
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

figure;
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


subplot(4,4,7)
a = [ABDr.NumberOfSaccadic]./[ABDr.NumberOfVestibular];
b = [ABDc.NumberOfSaccadic]./[ABDc.NumberOfVestibular];
c = [ABDIr.NumberOfSaccadic]./[ABDIr.NumberOfVestibular]; 
d = [ABDIc.NumberOfSaccadic]./[ABDIc.NumberOfVestibular];
    
plot(ones(1,size([ABDr.NumberOfRemaining],2)),a,'ko');
hold on;
plot(2*ones(1,size([ABDc.NumberOfRemaining],2)),b,'ko');
plot(3*ones(1,size([ABDIr.NumberOfRemaining],2)),c,'ko')
plot(4*ones(1,size([ABDIc.NumberOfRemaining],2)),d,'ko')

plot([1,2,3,4], [mean(a(~isinf(a))),mean(b(~isinf(b))), mean(c(~isinf(c))),mean(d(~isinf(d)))], '-ko', 'LineWidth',2);

axis square;
box off
set(gca, 'XTick',[1,2,3,4], 'YLim',[0,30],'XTickLabel',{'ABDr', 'ABDc', 'ABDIr', 'ABDIc'}, 'XTickLabelRotation',45);
ylabel('Saccade/Vestibular ratio');



%% plot positon of synapses on Z-brain space

% ABDr
for i = 1:numel(ABDr_CellIDs)
    ABDr(i).SaccadicDist = ABDr(i).PathLength(ABDr(i).isSaccadic)/max(Pvec_tree(ABDr(i).Tree{1}));
    ABDr(i).VestibularDist = ABDr(i).PathLength(ABDr(i).isVestibular)/max(Pvec_tree(ABDr(i).Tree{1}));
    ABDr(i).ContraDist = ABDr(i).PathLength(ABDr(i).isContra)/max(Pvec_tree(ABDr(i).Tree{1}));
    ABDr(i).IntegratorDist = ABDr(i).PathLength(ABDr(i).isIntegrator)/max(Pvec_tree(ABDr(i).Tree{1}));
    ABDr(i).EverythingElseDist = ABDr(i).PathLength(ABDr(i).isEverythingElse)/max(Pvec_tree(ABDr(i).Tree{1}));
end

% ABDc
for i = 1:numel(ABDc_CellIDs)
    if ~isempty(ABDc(i).Tree)
    ABDc(i).SaccadicDist  = ABDc(i).PathLength(ABDc(i).isSaccadic)/max(Pvec_tree(ABDc(i).Tree{1}));
    ABDc(i).VestibularDist  = ABDc(i).PathLength(ABDc(i).isVestibular)/max(Pvec_tree(ABDc(i).Tree{1}));
    ABDc(i).ContraDist  = ABDc(i).PathLength(ABDc(i).isContra)/max(Pvec_tree(ABDc(i).Tree{1}));
    ABDc(i).IntegratorDist  =  ABDc(i).PathLength(ABDc(i).isIntegrator)/max(Pvec_tree(ABDc(i).Tree{1}));
    ABDc(i).EverythingElseDist = ABDc(i).PathLength(ABDc(i).isEverythingElse)/max(Pvec_tree(ABDc(i).Tree{1}));
    end
end

% ABDIr

for i = 1:numel(ABDIr_CellIDs)
    ABDIr(i).SaccadicDist  = ABDIr(i).PathLength(ABDIr(i).isSaccadic)/max(Pvec_tree(ABDIr(i).Tree{1}));
    ABDIr(i).VestibularDist  = ABDIr(i).PathLength(ABDIr(i).isVestibular)/max(Pvec_tree(ABDIr(i).Tree{1}));
    ABDIr(i).ContraDist  = ABDIr(i).PathLength(ABDIr(i).isContra)/max(Pvec_tree(ABDIr(i).Tree{1}));
    ABDIr(i).IntegratorDist  = ABDIr(i).PathLength(ABDIr(i).isIntegrator)/max(Pvec_tree(ABDIr(i).Tree{1}));
    ABDIr(i).EverythingElseDist = ABDIr(i).PathLength(ABDIr(i).isEverythingElse)/max(Pvec_tree(ABDIr(i).Tree{1}));
end


% ABDIc

for i = 1:numel(ABDIc_CellIDs)
    ABDIc(i).SaccadicDist  = ABDIc(i).PathLength(ABDIc(i).isSaccadic)/max(Pvec_tree(ABDIc(i).Tree{1}));
    ABDIc(i).VestibularDist  = ABDIc(i).PathLength(ABDIc(i).isVestibular)/max(Pvec_tree(ABDIc(i).Tree{1}));
    ABDIc(i).ContraDist  = ABDIc(i).PathLength(ABDIc(i).isContra)/max(Pvec_tree(ABDIc(i).Tree{1}));
    ABDIc(i).IntegratorDist  = ABDIc(i).PathLength(ABDIc(i).isIntegrator)/max(Pvec_tree(ABDIc(i).Tree{1}));
    ABDIc(i).EverythingElseDist = ABDIc(i).PathLength(ABDIc(i).isEverythingElse)/max(Pvec_tree(ABDIc(i).Tree{1}));
end


getABDgradient(ABDc,ABDPutativeSaccadic.cellIDs,true)
%% Plot for each cell

for i = 1:numel(ABDr_CellIDs)
    subplot(4,4,i)
    histogram(ABDr(i).SaccadicDist,20,'FaceColor','none','EdgeColor',colors(1,:),'DisplayStyle','stairs','LineWidth',2);
    hold on;
    histogram(ABDr(i).VestibularDist,20,'FaceColor','none','EdgeColor',colors(2,:),'DisplayStyle','stairs','LineWidth',2);
    histogram(ABDr(i).ContraDist,20,'FaceColor','none','EdgeColor',colors(3,:),'DisplayStyle','stairs','LineWidth',2);
    histogram(ABDr(i).IntegratorDist,20,'FaceColor','none','EdgeColor',colors(4,:),'DisplayStyle','stairs','LineWidth',2);
    histogram(ABDr(i).EverythingElseDist,20,'FaceColor','none','EdgeColor',colors(5,:),'DisplayStyle','stairs','LineWidth',2);
    set(gca,'XLim',[0,1]);
end

%% plot distributions of everything
figure;

% saccadic pop
subplot(4,5,1)
h1 = histogram(vertcat(ABDr.SaccadicDist),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(1,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,1],'YLim',[0,0.15]);
box off;
subplot(4,5,6)
h1 = histogram(vertcat(ABDc.SaccadicDist),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(1,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,1],'YLim',[0,0.15]);
box off;
subplot(4,5,11)
h1 = histogram(vertcat(ABDIr.SaccadicDist),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(1,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,1],'YLim',[0,0.15]);
box off;
subplot(4,5,16)
h1 = histogram(vertcat(ABDIc.SaccadicDist),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(1,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,1],'YLim',[0,0.15]);
box off;

% Vestibular pop

subplot(4,5,2)
h1=histogram(vertcat(ABDr.VestibularDist),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(2,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,1],'YLim',[0,0.15]);
box off;
subplot(4,5,7)
h1=histogram(vertcat(ABDc.VestibularDist),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(2,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,1],'YLim',[0,0.15]);
box off;
subplot(4,5,12)
h1=histogram(vertcat(ABDIr.VestibularDist),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(2,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,1],'YLim',[0,0.15]);
box off;
subplot(4,5,17)
h1=histogram(vertcat(ABDIc.VestibularDist),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(2,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,1],'YLim',[0,0.15]);
box off;

% Integrator pop
subplot(4,5,3)
h1=histogram(vertcat(ABDr.IntegratorDist),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(4,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,1],'YLim',[0,0.15]);
box off;
subplot(4,5,8)
h1=histogram(vertcat(ABDc.IntegratorDist),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(4,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,1],'YLim',[0,0.15]);
box off;
subplot(4,5,13)
h1=histogram(vertcat(ABDIr.IntegratorDist),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(4,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,1],'YLim',[0,0.15]);
box off;
subplot(4,5,18)
h1=histogram(vertcat(ABDIc.IntegratorDist),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(4,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,1],'YLim',[0,0.15]);
box off;

% Contra pop
subplot(4,5,4)
h1=histogram(vertcat(ABDr.ContraDist),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(3,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,1],'YLim',[0,0.15]);
box off;
subplot(4,5,9)
h1=histogram(vertcat(ABDc.ContraDist),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(3,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,1],'YLim',[0,0.15]);
box off;
subplot(4,5,14)
h1=histogram(vertcat(ABDIr.ContraDist),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(3,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,1],'YLim',[0,0.15]);
box off;
subplot(4,5,19)
h1=histogram(vertcat(ABDIc.ContraDist),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(3,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,1],'YLim',[0,0.15]);
box off;

% Remaining pop
subplot(4,5,5)
h1=histogram(vertcat(ABDr.EverythingElseDist),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(5,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,1],'YLim',[0,0.15]);
box off;
subplot(4,5,10)
h1=histogram(vertcat(ABDc.EverythingElseDist),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(5,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,1],'YLim',[0,0.15]);
box off;
subplot(4,5,15)
h1=histogram(vertcat(ABDIr.EverythingElseDist),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(5,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,1],'YLim',[0,0.15]);
box off;
subplot(4,5,20)
h1=histogram(vertcat(ABDIc.EverythingElseDist),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = colors(5,:);
h1.LineWidth = 2;
set(gca,'XLim',[0,1],'YLim',[0,0.15]);
box off;

%% Cumulative Plots

figure;

subplot(4,4,1);

colormap = cbrewer('qual','Set1',5);

histogram([ABDPutativeSaccadic.ABDcpathLength;ABDPutativeSaccadic.ABDrpathLength],...
    20,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',colormap(1,:)); % Saccadic
hold on;
histogram([vertcat(ABDr.SaccadicDist);vertcat(ABDc.SaccadicDist);vertcat(ABDr.IntegratorDist);vertcat(ABDc.IntegratorDist)],...
    20,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',colormap(2,:)) ; % r45,r67 combined Integrator
histogram([vertcat(ABDr.VestibularDist);vertcat(ABDc.VestibularDist)],...
    20,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',colormap(3,:));

legend({'S','Int','V','Contra'},'Location','bestoutside','Box','OFF');

set(gca,'XLim',[0,1],'YLim',[0,1]);
xlabel('Norm. pathlength');
ylabel('Cumulative count');
box off;
offsetAxes
axis square;

subplot(4,4,2);

histogram([ABDPutativeSaccadic.ABDcpathLength;ABDPutativeSaccadic.ABDrpathLength],...
    20,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',colormap(1,:)); % Saccadic
hold on;
histogram([vertcat(ABDr.SaccadicDist);vertcat(ABDc.SaccadicDist)],...
    20,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',colormap(2,:)) ; % r45,r67 combined Integrator

histogram([vertcat(ABDr.IntegratorDist);vertcat(ABDc.IntegratorDist)],...
    20,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',colormap(3,:)) ; % r45,r67 combined Integrator

histogram([vertcat(ABDr.VestibularDist);vertcat(ABDc.VestibularDist)],...
    20,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',colormap(4,:));
histogram([vertcat(ABDr.ContraDist);vertcat(ABDc.ContraDist)],...
    20,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',colormap(5,:));

legend({'S','r23','r78','V','Contra'},'Location','bestoutside','Box','OFF');

set(gca,'XLim',[0,1],'YLim',[0,1]);
xlabel('Norm. pathlength');
ylabel('Cumulative count');
box off;
offsetAxes
axis square;


subplot(4,4,3);

histogram([ABDPutativeSaccadic.ABDcpathLength;ABDPutativeSaccadic.ABDrpathLength],...
    20,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',colormap(1,:)); % Saccadic
hold on;
histogram([vertcat(ABDr.SaccadicDist);vertcat(ABDc.SaccadicDist);vertcat(ABDr.IntegratorDist);vertcat(ABDc.IntegratorDist)],...
    20,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',colormap(2,:)) ; % r45,r67 combined Integrator
histogram([vertcat(ABDr.VestibularDist);vertcat(ABDc.VestibularDist)],...
    20,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',colormap(3,:));
histogram([vertcat(ABDr.ContraDist);vertcat(ABDc.ContraDist)],...
    20,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',colormap(4,:));

legend({'S','Int','V','Contra'},'Location','bestoutside','Box','OFF');

set(gca,'XLim',[0,1],'YLim',[0,1]);
xlabel('Norm. pathlength');
ylabel('Cumulative count');
box off;
offsetAxes
axis square;

subplot(4,4,5);

histogram([ABDiPutativeSaccadic.ABDIcpathLength;ABDiPutativeSaccadic.ABDIrpathLength],...
    20,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',colormap(1,:)); % Saccadic
hold on;
histogram([vertcat(ABDIr.SaccadicDist);vertcat(ABDIc.SaccadicDist);vertcat(ABDIr.IntegratorDist);vertcat(ABDIc.IntegratorDist)],...
    20,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',colormap(2,:)) ; % r45,r67 combined Integrator
histogram([vertcat(ABDIr.VestibularDist);vertcat(ABDIc.VestibularDist)],...
    20,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',colormap(3,:));
% histogram([vertcat(ABDIr.ContraDist);vertcat(ABDIc.ContraDist)],...
%     20,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2);

legend({'S','Int','V','Contra'},'Location','bestoutside','Box','OFF');

set(gca,'XLim',[0,1],'YLim',[0,1]);
xlabel('Norm. pathlength');
ylabel('Cumulative count');
box off;
offsetAxes
axis square;

subplot(4,4,6);

%colormap = colorcet('R3','N',6);
histogram([ABDiPutativeSaccadic.ABDIcpathLength;ABDiPutativeSaccadic.ABDIrpathLength],...
    20,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',colormap(1,:)); % Saccadic
hold on;
histogram([vertcat(ABDIr.SaccadicDist);vertcat(ABDIc.SaccadicDist)],...
    20,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',colormap(2,:)) ; % r45,r67 combined Integrator

histogram([vertcat(ABDIr.IntegratorDist);vertcat(ABDIc.IntegratorDist)],...
    20,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',colormap(3,:)) ; % r45,r67 combined Integrator

histogram([vertcat(ABDIr.VestibularDist);vertcat(ABDIc.VestibularDist)],...
    20,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',colormap(4,:));
histogram([vertcat(ABDIr.ContraDist);vertcat(ABDIc.ContraDist)],...
    20,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',colormap(5,:));

legend({'S','r23','r78','V','Contra'},'Location','bestoutside','Box','OFF');

set(gca,'XLim',[0,1],'YLim',[0,1]);
xlabel('Norm. pathlength');
ylabel('Cumulative count');
box off;
offsetAxes
axis square;


subplot(4,4,7);

%colormap = colorcet('R3','N',5);
histogram([ABDiPutativeSaccadic.ABDIcpathLength;ABDiPutativeSaccadic.ABDIrpathLength],...
    20,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',colormap(1,:)); % Saccadic
hold on;
histogram([vertcat(ABDIr.SaccadicDist);vertcat(ABDIc.SaccadicDist);vertcat(ABDIr.IntegratorDist);vertcat(ABDIc.IntegratorDist)],...
    20,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',colormap(2,:)) ; % r45,r67 combined Integrator
histogram([vertcat(ABDIr.VestibularDist);vertcat(ABDIc.VestibularDist)],...
    20,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',colormap(3,:));
histogram([vertcat(ABDIr.ContraDist);vertcat(ABDIc.ContraDist)],...
    20,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',colormap(4,:));

legend({'S','Int','V','Contra'},'Location','bestoutside','Box','OFF');

set(gca,'XLim',[0,1],'YLim',[0,1]);
xlabel('Norm. pathlength');
ylabel('Cumulative count');
box off;
offsetAxes
axis square;



%% distributions of each class onto lumped Motor classes.

figure; 
% Saccadic pop
subplot(4,4,1)
h1=histogram([vertcat(ABDr.SaccadicDist);vertcat(ABDc.SaccadicDist)],20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = 'm';
h1.LineWidth = 2;
set(gca,'XLim',[0,1],'YLim',[0,0.15]);
box off;
hold on;

h1=histogram([vertcat(ABDIr.SaccadicDist);vertcat(ABDIc.SaccadicDist)],20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = 'g';
h1.LineWidth = 2;
set(gca,'XLim',[0,1],'YLim',[0,0.15]);
box off;
hold on;
axis square
xlabel('Pathlength (um)');
title('Saccadic axons');

% Vestibular pop
subplot(4,4,2)
h1=histogram([vertcat(ABDr.VestibularDist);vertcat(ABDc.VestibularDist)],20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = 'm';
h1.LineWidth = 2;
set(gca,'XLim',[0,1],'YLim',[0,0.15]);
box off;
hold on;

h1=histogram([vertcat(ABDIr.VestibularDist);vertcat(ABDIc.VestibularDist)],20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = 'g';
h1.LineWidth = 2;
set(gca,'XLim',[0,1],'YLim',[0,0.15]);
box off;
hold on;
axis square
title('Vestibular axons');

% Integrator pop
subplot(4,4,3)
h1=histogram([vertcat(ABDr.IntegratorDist);vertcat(ABDc.IntegratorDist)],20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = 'm';
h1.LineWidth = 2;
set(gca,'XLim',[0,1],'YLim',[0,0.15]);
box off;
hold on;

h1=histogram([vertcat(ABDIr.IntegratorDist);vertcat(ABDIc.IntegratorDist)],20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = 'g';
h1.LineWidth = 2;
set(gca,'XLim',[0,1],'YLim',[0,0.15]);
box off;
hold on;
axis square
title('Integrator axons');

%Contra pop
subplot(4,4,4)
h1=histogram([vertcat(ABDr.ContraDist);vertcat(ABDc.ContraDist)],20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = 'm';
h1.LineWidth = 2;
set(gca,'XLim',[0,1],'YLim',[0,0.15]);
box off;
hold on;

h1=histogram([vertcat(ABDIr.ContraDist);vertcat(ABDIc.ContraDist)],20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = 'g';
h1.LineWidth = 2;
set(gca,'XLim',[0,1],'YLim',[0,0.15]);
box off;
hold on;
axis square
title('Contra axons');

% Remaining pop
subplot(4,4,5)
h1=histogram([vertcat(ABDr.EverythingElseDist);vertcat(ABDc.EverythingElseDist)],20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = 'm';
h1.LineWidth = 2;
set(gca,'XLim',[0,1],'YLim',[0,0.15]);
box off;
hold on;

h1=histogram([vertcat(ABDIr.EverythingElseDist);vertcat(ABDIc.EverythingElseDist)],20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = 'g';
h1.LineWidth = 2;
set(gca,'XLim',[0,1],'YLim',[0,0.15]);
box off;
hold on;
axis square
title('Remaining axons');

legend({'ABD','ABDi'},'Location','bestoutside');

%% Split the contributions from each of the axons classes

load('leadSynapseDiff.mat','leadSynapseDiff')
load('leadMotorNeuronDiff.mat','leadMotorNeuronDiff');
load('leadDiffAxons.mat','leadDiffAxons');

load('lagSynapseDiff.mat','lagSynapseDiff')
load('leadMotorNeuronDiff.mat','leadMotorNeuronDiff');
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
    [~,ia,~] = intersect(ABDr(i).Saccadic, leadDiffAxons.Saccadic);
    for j = 1:length(ia)
    tempdist = ABDr(i).SaccadicDist(ia(j));
    %tempdist = PathLengthToCoordinate(ABDr(i).PreSynCoordsTransformed(ia(j),:),ABDr(i).Tree{1});
    temp = [temp; tempdist];
    end
    clear ia;
    ABDr(i).LeadSaccadicDist = temp;
    temp = [];
    
    [~,ia,~] = intersect(ABDr(i).Saccadic, [lagDiffAxons.Saccadic]);
    for j = 1:length(ia)
        tempdist = ABDr(i).SaccadicDist(ia(j));
   % tempdist = PathLengthToCoordinate(ABDr(i).PreSynCoordsTransformed(ia(j),:),ABDr(i).Tree{1});
    temp = [temp; tempdist];
    end
    clear ia;
    ABDr(i).LagSaccadicDist = temp;
    temp =[];
end

temp = [];
for i = 1:numel(ABDc_CellIDs)
    [~,ia,~] = intersect(ABDc(i).Saccadic, [leadDiffAxons.Saccadic]);
    for j = 1:length(ia)
    tempdist = ABDc(i).SaccadicDist(ia(j));
    temp = [temp; tempdist];
    end
    clear ia;
    ABDc(i).LeadSaccadicDist = temp;
    temp = [];
    
    [~,ia,~] = intersect(ABDc(i).Saccadic, [lagDiffAxons.Saccadic]);
    for j = 1:length(ia)
    tempdist = ABDc(i).SaccadicDist(ia(j));
    temp = [temp; tempdist];
    end
    clear ia;
    ABDc(i).LagSaccadicDist = temp;
    temp =[];
end
    

subplot(4,4,2)
h1=histogram([vertcat(ABDr.LeadSaccadicDist);vertcat(ABDc.LeadSaccadicDist)],20,'FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = leadColor;
h1.LineWidth = 2;
%set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;
h1=histogram([vertcat(ABDr.LagSaccadicDist);vertcat(ABDc.LagSaccadicDist)],20,'FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = lagColor;
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
    [~,ia,~] = intersect(ABDIr(i).Saccadic, [leadDiffAxons.Saccadic]);
    for j = 1:length(ia)
    tempdist = ABDIr(i).SaccadicDist(ia(j));
    temp = [temp; tempdist];
    end
    clear ia;
    ABDIr(i).LeadSaccadicDist = temp;
    temp = [];
    
    [~,ia,~] = intersect(ABDIr(i).Saccadic, [lagDiffAxons.Saccadic]);
    for j = 1:length(ia)
    tempdist = ABDIr(i).SaccadicDist(ia(j));
    temp = [temp; tempdist];
    end
    clear ia;
    ABDIr(i).LagSaccadicDist = temp;
    temp =[];
end

temp = [];
for i = 1:numel(ABDIc_CellIDs)
    [~,ia,~] = intersect(ABDIc(i).Saccadic, [leadDiffAxons.Saccadic]);
    for j = 1:length(ia)
    tempdist = ABDIc(i).SaccadicDist(ia(j));
    temp = [temp; tempdist];
    end
    clear ia;
    ABDIc(i).LeadSaccadicDist = temp;
    temp = [];
    
    [~,ia,~] = intersect(ABDIc(i).Saccadic, [lagDiffAxons.Saccadic]);
    for j = 1:length(ia)
    tempdist = ABDIc(i).SaccadicDist(ia(j));
    temp = [temp; tempdist];
    end
    clear ia;
    ABDIc(i).LagSaccadicDist = temp;
    temp =[];
end
    

subplot(4,4,3)
h1=histogram([vertcat(ABDIr.LeadSaccadicDist);vertcat(ABDIc.LeadSaccadicDist)],20,'FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = leadColor;
h1.LineWidth = 2;
%set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;
h1=histogram([vertcat(ABDIr.LagSaccadicDist);vertcat(ABDIc.LagSaccadicDist)],20,'FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = lagColor;
h1.LineWidth = 2;
%set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;
axis square
xlabel('Pathlength (um)');
title('ABDi Saccadic axons');

subplot(4,4,4)
SaccadicAxons = [vertcat(ABDr.Saccadic);vertcat(ABDc.Saccadic);vertcat(ABDIr.Saccadic);vertcat(ABDIc.Saccadic)];
uniqueSaccadicAxons = unique(SaccadicAxons);
saccadicMotorSynapses = isMotor(uniqueSaccadicAxons,df);
saccadicMotorNeuronTargets = isPostSynapseMotor(uniqueSaccadicAxons,df);

scatter([sum(saccadicMotorSynapses(:,2:3),2)-sum(saccadicMotorSynapses(:,4:5),2)], ...
    [sum(saccadicMotorNeuronTargets(:,1:2),2)-sum(saccadicMotorNeuronTargets(:,3:4),2)],20,'o',...
    'MarkerFaceColor',colors(1,:),'MarkerEdgeColor','none');
line([0,0],[-20,20],'linestyle',':','color','k');
line([-200,200],[0,0],'linestyle',':','color','k');
ylabel('#(ABD - ABDi) neurons');
xlabel('#(ABD - ABDi) synapses')
title('Saccadic Axons');
text(100,-18,sprintf('n = %d',size(uniqueSaccadicAxons,1)));
box off;
axis square;


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
title('Vestibular axons');

temp = [];
for i = 1:numel(ABDr_CellIDs)
    [~,ia,~] = intersect(ABDr(i).Vestibular, [leadDiffAxons.Vestibular]);
    for j = 1:length(ia)
    tempdist = ABDr(i).VestibularDist(ia(j));
    temp = [temp; tempdist];
    end
    clear ia;
    ABDr(i).LeadVestibularDist = temp;
    temp = [];
    
    [~,ia,~] = intersect(ABDr(i).Vestibular, [lagDiffAxons.Vestibular]);
    for j = 1:length(ia)
    tempdist = ABDr(i).VestibularDist(ia(j));
    temp = [temp; tempdist];
    end
    clear ia;
    ABDr(i).LagVestibularDist = temp;
    temp =[];
end

temp = [];
for i = 1:numel(ABDc_CellIDs)
    [~,ia,~] = intersect(ABDc(i).Vestibular, [leadDiffAxons.Vestibular]);
    for j = 1:length(ia)
    tempdist = ABDc(i).VestibularDist(ia(j));
    temp = [temp; tempdist];
    end
    clear ia;
    ABDc(i).LeadVestibularDist = temp;
    temp = [];
    
    [~,ia,~] = intersect(ABDc(i).Vestibular, [lagDiffAxons.Vestibular]);
    for j = 1:length(ia)
    tempdist = ABDc(i).VestibularDist(ia(j));
    temp = [temp; tempdist];
    end
    clear ia;
    ABDc(i).LagVestibularDist = temp;
    temp =[];
end
    

subplot(4,4,6)
h1=histogram([vertcat(ABDr.LeadVestibularDist);vertcat(ABDc.LeadVestibularDist)],20,'FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = leadColor;
h1.LineWidth = 2;
%set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;
h1=histogram([vertcat(ABDr.LagVestibularDist);vertcat(ABDc.LagVestibularDist)],20,'FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = lagColor;
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
    [~,ia,~] = intersect(ABDIr(i).Vestibular, [leadDiffAxons.Vestibular]);
    for j = 1:length(ia)
    tempdist = ABDIr(i).VestibularDist(ia(j));
    temp = [temp; tempdist];
    end
    clear ia;
    ABDIr(i).LeadVestibularDist = temp;
    temp = [];
    
    [~,ia,~] = intersect(ABDIr(i).Vestibular, [lagDiffAxons.Vestibular]);
    for j = 1:length(ia)
    tempdist = ABDIr(i).VestibularDist(ia(j));
    temp = [temp; tempdist];
    end
    clear ia;
    ABDIr(i).LagVestibularDist = temp;
    temp =[];
end

temp = [];
for i = 1:numel(ABDIc_CellIDs)
    [~,ia,~] = intersect(ABDIc(i).Vestibular, [leadDiffAxons.Vestibular]);
    for j = 1:length(ia)
    tempdist = ABDIc(i).VestibularDist(ia(j));
    temp = [temp; tempdist];
    end
    clear ia;
    ABDIc(i).LeadVestibularDist = temp;
    temp = [];
    
    [~,ia,~] = intersect(ABDIc(i).Vestibular, [lagDiffAxons.Vestibular]);
    for j = 1:length(ia)
    tempdist = ABDIc(i).VestibularDist(ia(j));
    temp = [temp; tempdist];
    end
    clear ia;
    ABDIc(i).LagVestibularDist = temp;
    temp =[];
end
    

subplot(4,4,7)
h1=histogram([vertcat(ABDIr.LeadVestibularDist);vertcat(ABDIc.LeadVestibularDist)],20,'FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = leadColor;
h1.LineWidth = 2;
%set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;
h1=histogram([vertcat(ABDIr.LagVestibularDist);vertcat(ABDIc.LagVestibularDist)],20,'FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = lagColor;
h1.LineWidth = 2;
%set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;
axis square
xlabel('Pathlength (um)');
title('ABDi Vestibular axons');

subplot(4,4,8)
VestibularAxons = [vertcat(ABDr.Vestibular);vertcat(ABDc.Vestibular);vertcat(ABDIr.Vestibular);vertcat(ABDIc.Vestibular)];
uniqueVestibularAxons = unique(VestibularAxons);
vestibularMotorSynapses = isMotor(uniqueVestibularAxons,df);
vestibularMotorNeuronTargets = isPostSynapseMotor(uniqueVestibularAxons,df);

scatter([sum(vestibularMotorSynapses(:,2:3),2)-sum(vestibularMotorSynapses(:,4:5),2)], ...
    [sum(vestibularMotorNeuronTargets(:,1:2),2)-sum(vestibularMotorNeuronTargets(:,3:4),2)],20,'o',...
    'MarkerFaceColor',colors(2,:),'MarkerEdgeColor','none');
line([0,0],[-20,20],'linestyle',':','color','k');
line([-200,200],[0,0],'linestyle',':','color','k');
ylabel('#(ABD - ABDi) neurons');
xlabel('#(ABD - ABDi) synapses');
text(100,-18,sprintf('n = %d',size(uniqueVestibularAxons,1)));
title('Vestibular Axons');
box off;
axis square;

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
    [~,ia,~] = intersect(ABDr(i).Integrator, [leadDiffAxons.Integrator]);
    for j = 1:length(ia)
    tempdist = ABDr(i).IntegratorDist(ia(j));
    temp = [temp; tempdist];
    end
    clear ia;
    ABDr(i).LeadIntegratorDist = temp;
    temp = [];
    
    [~,ia,~] = intersect(ABDr(i).Integrator, [lagDiffAxons.Integrator]);
    for j = 1:length(ia)
    tempdist = ABDr(i).IntegratorDist(ia(j));
    temp = [temp; tempdist];
    end
    clear ia;
    ABDr(i).LagIntegratorDist = temp;
    temp =[];
end

temp = [];
for i = 1:numel(ABDc_CellIDs)
    [~,ia,~] = intersect(ABDc(i).Integrator, [leadDiffAxons.Integrator]);
    for j = 1:length(ia)
    tempdist = ABDc(i).IntegratorDist(ia(j));
    temp = [temp; tempdist];
    end
    clear ia;
    ABDc(i).LeadIntegratorDist = temp;
    temp = [];
    
    [~,ia,~] = intersect(ABDc(i).Integrator, [lagDiffAxons.Integrator]);
    for j = 1:length(ia)
    tempdist = ABDc(i).IntegratorDist(ia(j));
    temp = [temp; tempdist];
    end
    clear ia;
    ABDc(i).LagIntegratorDist = temp;
    temp =[];
end
    

subplot(4,4,10)
h1=histogram([vertcat(ABDr.LeadIntegratorDist);vertcat(ABDc.LeadIntegratorDist)],20,'FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = leadColor;
h1.LineWidth = 2;
%set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;
h1=histogram([vertcat(ABDr.LagIntegratorDist);vertcat(ABDc.LagIntegratorDist)],20,'FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = lagColor;
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
    [~,ia,~] = intersect(ABDIr(i).Integrator, [leadDiffAxons.Integrator]);
    for j = 1:length(ia)  
    tempdist = ABDIr(i).IntegratorDist(ia(j));
    temp = [temp; tempdist];
    end
    clear ia;
    ABDIr(i).LeadIntegratorDist = temp;
    temp = [];
    
    [~,ia,~] = intersect(ABDIr(i).Integrator, [lagDiffAxons.Integrator]);
    for j = 1:length(ia)
    tempdist = ABDIr(i).IntegratorDist(ia(j));
    temp = [temp; tempdist];
    end
    clear ia;
    ABDIr(i).LagIntegratorDist = temp;
    temp =[];
end

temp = [];
for i = 1:numel(ABDIc_CellIDs)
    [~,ia,~] = intersect(ABDIc(i).Integrator, [leadDiffAxons.Integrator]);
    for j = 1:length(ia)
    tempdist = ABDIc(i).IntegratorDist(ia(j));
    temp = [temp; tempdist];
    end
    clear ia;
    ABDIc(i).LeadIntegratorDist = temp;
    temp = [];
    
    [~,ia,~] = intersect(ABDIc(i).Integrator, [lagDiffAxons.Integrator]);
    for j = 1:length(ia)
    tempdist = ABDIc(i).IntegratorDist(ia(j));
    temp = [temp; tempdist];
    end
    clear ia;
    ABDIc(i).LagIntegratorDist = temp;
    temp =[];
end
    

subplot(4,4,11)
h1=histogram([vertcat(ABDIr.LeadIntegratorDist);vertcat(ABDIc.LeadIntegratorDist)],20,'FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = leadColor;
h1.LineWidth = 2;
%set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;
h1=histogram([vertcat(ABDIr.LagIntegratorDist);vertcat(ABDIc.LagIntegratorDist)],20,'FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = lagColor;
h1.LineWidth = 2;
%set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;
axis square
xlabel('Pathlength (um)');
title('ABDi Integrator axons');

subplot(4,4,12)
IntegratorAxons = [vertcat(ABDr.Integrator);vertcat(ABDc.Integrator);vertcat(ABDIr.Integrator);vertcat(ABDIc.Integrator)];
uniqueIntegratorAxons = unique(IntegratorAxons);
integratorMotorSynapses = isMotor(uniqueIntegratorAxons,df);
integratorMotorNeuronTargets = isPostSynapseMotor(uniqueIntegratorAxons,df);

scatter([sum(integratorMotorSynapses(:,2:3),2)-sum(integratorMotorSynapses(:,4:5),2)], ...
    [sum(integratorMotorNeuronTargets(:,1:2),2)-sum(integratorMotorNeuronTargets(:,3:4),2)],20,'o',...
    'MarkerFaceColor',colors(3,:),'MarkerEdgeColor','none');
line([0,0],[-20,20],'linestyle',':','color','k');
line([-200,200],[0,0],'linestyle',':','color','k');
ylabel('#(ABD - ABDi) neurons');
xlabel('#(ABD - ABDi) synapses')
text(100,-18,sprintf('n = %d',size(uniqueIntegratorAxons,1)));
title('Integrator Axons');
box off;
axis square;

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
    [~,ia,~] = intersect(ABDr(i).Contra, [leadDiffAxons.Contra]);
    for j = 1:length(ia)
    tempdist = ABDr(i).ContraDist(ia(j));
    temp = [temp; tempdist];
    end
    clear ia;
    ABDr(i).LeadContraDist = temp;
    temp = [];
    
    [~,ia,~] = intersect(ABDr(i).Contra, [lagDiffAxons.Contra]);
    for j = 1:length(ia)
    tempdist = ABDr(i).ContraDist(ia(j));
    temp = [temp; tempdist];
    end
    clear ia;
    ABDr(i).LagContraDist = temp;
    temp =[];
end

temp = [];
for i = 1:numel(ABDc_CellIDs)
    [~,ia,~] = intersect(ABDc(i).Contra, [leadDiffAxons.Contra]);
    for j = 1:length(ia)
    tempdist = ABDc(i).ContraDist(ia(j));
    temp = [temp; tempdist];
    end
    clear ia;
    ABDc(i).LeadContraDist = temp;
    temp = [];
    
    [~,ia,~] = intersect(ABDc(i).Contra, [lagDiffAxons.Contra]);
    for j = 1:length(ia)
    tempdist = ABDc(i).ContraDist(ia(j));
    temp = [temp; tempdist];
    end
    clear ia;
    ABDc(i).LagContraDist = temp;
    temp =[];
end
    

subplot(4,4,14)
h1=histogram([vertcat(ABDr.LeadContraDist);vertcat(ABDc.LeadContraDist)],20,'FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = leadColor;
h1.LineWidth = 2;
%set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;
h1=histogram([vertcat(ABDr.LagContraDist);vertcat(ABDc.LagContraDist)],20,'FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = lagColor;
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
    [~,ia,~] = intersect(ABDIr(i).Contra, [leadDiffAxons.Contra]);
    for j = 1:length(ia)
    tempdist = ABDIr(i).ContraDist(ia(j));
    temp = [temp; tempdist];
    end
    clear ia;
    ABDIr(i).LeadContraDist = temp;
    temp = [];
    
    [~,ia,~] = intersect(ABDIr(i).Contra, [lagDiffAxons.Contra]);
    for j = 1:length(ia)
    tempdist = ABDIr(i).ContraDist(ia(j));
    temp = [temp; tempdist];
    end
    clear ia;
    ABDIr(i).LagContraDist = temp;
    temp =[];
end

temp = [];
for i = 1:numel(ABDIc_CellIDs)
    [~,ia,~] = intersect(ABDIc(i).Contra, [leadDiffAxons.Contra]);
    for j = 1:length(ia)
    tempdist = ABDIc(i).ContraDist(ia(j));
    temp = [temp; tempdist];
    end
    clear ia;
    ABDIc(i).LeadContraDist = temp;
    temp = [];
    
    [~,ia,~] = intersect(ABDIc(i).Contra, [lagDiffAxons.Contra]);
    for j = 1:length(ia)
    tempdist = ABDIc(i).ContraDist(ia(j));
    temp = [temp; tempdist];
    end
    clear ia;
    ABDIc(i).LagContraDist = temp;
    temp =[];
end
    

subplot(4,4,15)
h1=histogram([vertcat(ABDIr.LeadContraDist);vertcat(ABDIc.LeadContraDist)],20,'FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = leadColor;
h1.LineWidth = 2;
%set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;
h1=histogram([vertcat(ABDIr.LagContraDist);vertcat(ABDIc.LagContraDist)],20,'FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = lagColor;
h1.LineWidth = 2;
%set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;
axis square
xlabel('Pathlength (um)');
title('ABDi Contra axons');

subplot(4,4,16)
ContraAxons = [vertcat(ABDr.Contra);vertcat(ABDc.Contra);vertcat(ABDIr.Contra);vertcat(ABDIc.Contra)];
uniqueContraAxons = unique(ContraAxons);
contraMotorSynapses = isMotor(uniqueContraAxons,df);
contraMotorNeuronTargets = isPostSynapseMotor(uniqueContraAxons,df);

scatter([sum(contraMotorSynapses(:,2:3),2)-sum(contraMotorSynapses(:,4:5),2)], ...
    [sum(contraMotorNeuronTargets(:,1:2),2)-sum(contraMotorNeuronTargets(:,3:4),2)],20,'o',...
    'MarkerFaceColor',colors(4,:),'MarkerEdgeColor','none');
line([0,0],[-20,20],'linestyle',':','color','k');
line([-200,200],[0,0],'linestyle',':','color','k');
ylabel('#(ABD - ABDi) neurons');
xlabel('#(ABD - ABDi) synapses')
title('Contra Axons');
text(100,-18,sprintf('n = %d',size(uniqueContraAxons,1)));
box off;
axis square;

%% plot the bimodal contra distribution

figure;
Contra_ABD_ABDi = [contraMotorNeuronTargets(:,1), [sum(contraMotorNeuronTargets(:,2:3),2)-sum(contraMotorNeuronTargets(:,4:5),2)]];
subplot(4,4,1);
histogram(Contra_ABD_ABDi(:,2));
axis square;
box off;

% set threshold to 10

Contra_ABD_heavy = Contra_ABD_ABDi(find(Contra_ABD_ABDi(:,2)>10));
Contra_ABDi_heavy = Contra_ABD_ABDi(find(Contra_ABD_ABDi(:,2)<-10));

Contra_ABD_rem = Contra_ABD_ABDi(find(Contra_ABD_ABDi(:,2)>0 & Contra_ABD_ABDi(:,2)<10));
Contra_ABDi_rem = Contra_ABD_ABDi(find(Contra_ABD_ABDi(:,2)<0 & Contra_ABD_ABDi(:,2)>-10));



A = colorcet('D1','N',256,'reverse',1);
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(1,2,1);
transform_swc_AV(Contra_ABD_heavy,A(1,:),[],false,false);
transform_swc_AV(Contra_ABDi_heavy,A(256,:),[],false,false);

subplot(1,2,2)
transform_swc_AV(Contra_ABD_rem,A(1,:),[],true,false);
transform_swc_AV(Contra_ABDi_rem,A(256,:),[],false,false);



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
    [~,ia,~] = intersect(ABDr(i).EverythingElse, [leadDiffAxons.EverythingElse]);
    for j = 1:length(ia)
    tempdist = ABDr(i).EverythingElseDist(ia(j));
    temp = [temp; tempdist];
    end
    clear ia;
    ABDr(i).LeadEverythingElseDist = temp;
    temp = [];
    
    [~,ia,~] = intersect(ABDr(i).EverythingElse, [lagDiffAxons.EverythingElse]);
    for j = 1:length(ia)
    tempdist = ABDr(i).EverythingElseDist(ia(j));
    temp = [temp; tempdist];
    end
    clear ia;
    ABDr(i).LagEverythingElseDist = temp;
    temp =[];
end

temp = [];
for i = 1:numel(ABDc_CellIDs)
    [~,ia,~] = intersect(ABDc(i).EverythingElse, [leadDiffAxons.EverythingElse]);
    for j = 1:length(ia)
    tempdist = ABDc(i).EverythingElseDist(ia(j));
    temp = [temp; tempdist];
    end
    clear ia;
    ABDc(i).LeadEverythingElseDist = temp;
    temp = [];
    
    [~,ia,~] = intersect(ABDc(i).EverythingElse, [lagDiffAxons.EverythingElse]);
    for j = 1:length(ia)
    tempdist = ABDc(i).EverythingElseDist(ia(j));
    temp = [temp; tempdist];
    end
    clear ia;
    ABDc(i).LagEverythingElseDist = temp;
    temp =[];
end
    

subplot(4,4,2)
h1=histogram([vertcat(ABDr.LeadEverythingElseDist);vertcat(ABDc.LeadEverythingElseDist)],20,'FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = leadColor;
h1.LineWidth = 2;
%set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;
h1=histogram([vertcat(ABDr.LagEverythingElseDist);vertcat(ABDc.LagEverythingElseDist)],20,'FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = lagColor;
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
    [~,ia,~] = intersect(ABDIr(i).EverythingElse, [leadDiffAxons.EverythingElse]);
    for j = 1:length(ia)
    tempdist = ABDIr(i).EverythingElseDist(ia(j));
    temp = [temp; tempdist];
    end
    clear ia;
    ABDIr(i).LeadEverythingElseDist = temp;
    temp = [];
    
    [~,ia,~] = intersect(ABDIr(i).EverythingElse, [lagDiffAxons.EverythingElse]);
    for j = 1:length(ia)
    tempdist = ABDIr(i).EverythingElseDist(ia(j));
    temp = [temp; tempdist];
    end
    clear ia;
    ABDIr(i).LagEverythingElseDist = temp;
    temp =[];
end

temp = [];
for i = 1:numel(ABDIc_CellIDs)
    [~,ia,~] = intersect(ABDIc(i).EverythingElse, [leadDiffAxons.EverythingElse]);
    for j = 1:length(ia)
    tempdist = ABDIc(i).EverythingElseDist(ia(j));
    temp = [temp; tempdist];
    end
    clear ia;
    ABDIc(i).LeadEverythingElseDist = temp;
    temp = [];
    
    [~,ia,~] = intersect(ABDIc(i).EverythingElse, [lagDiffAxons.EverythingElse]);
    for j = 1:length(ia)
    tempdist = ABDIc(i).EverythingElseDist(ia(j));
    temp = [temp; tempdist];
    end
    clear ia;
    ABDIc(i).LagEverythingElseDist = temp;
    temp =[];
end
    

subplot(4,4,3)
h1=histogram([vertcat(ABDIr.LeadEverythingElseDist);vertcat(ABDIc.LeadEverythingElseDist)],20,'FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = leadColor;
h1.LineWidth = 2;
%set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;
h1=histogram([vertcat(ABDIr.LagEverythingElseDist);vertcat(ABDIc.LagEverythingElseDist)],20,'FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = lagColor;
h1.LineWidth = 2;
%set(gca,'XLim',[0,100],'YLim',[0,0.15]);
box off;
hold on;
axis square
xlabel('Pathlength (um)');
title('ABDi EverythingElse axons');

subplot(4,4,4)
EverythingElseAxons = [vertcat(ABDr.EverythingElse);vertcat(ABDc.EverythingElse);vertcat(ABDIr.EverythingElse);vertcat(ABDIc.EverythingElse)];
uniqueEverythingElseAxons = unique(EverythingElseAxons);
uniqueEverythingElseAxons = uniqueEverythingElseAxons(uniqueEverythingElseAxons<1e5);
EverythingElseMotorSynapses = isMotor(uniqueEverythingElseAxons,df);
EverythingElseAxonsMotorNeuronTargets = isPostSynapseMotor(uniqueEverythingElseAxons,df);

scatter([sum(EverythingElseMotorSynapses(:,2:3),2)-sum(EverythingElseMotorSynapses(:,4:5),2)], ...
    [sum(EverythingElseAxonsMotorNeuronTargets(:,1:2),2)-sum(EverythingElseAxonsMotorNeuronTargets(:,3:4),2)],20,'o',...
    'MarkerFaceColor',colors(5,:),'MarkerEdgeColor','none');
line([0,0],[-20,20],'linestyle',':','color','k');
line([-200,200],[0,0],'linestyle',':','color','k');
ylabel('#(ABD - ABDi) neurons');
xlabel('#(ABD - ABDi) synapses')
title('Remaining Axons');
text(100,-18,sprintf('n = %d',size(uniqueEverythingElseAxons,1)));
box off;
axis square;
%% plot the skeleting of all the different classes the project to the lead/lag pop
A = colorcet('D1','N',256,'reverse',1);

% Saccadic pop
figure('units','normalized','outerposition',[0 0 1 1]);
temp = [sum(saccadicMotorNeuronTargets(:,2:3),2)-sum(saccadicMotorNeuronTargets(:,4:5),2)];
SaccadicOntoABD = uniqueSaccadicAxons(temp>0);
SaccadicOntoABDi = uniqueSaccadicAxons(temp<0);

subplot(1,2,1)
transform_swc_AV([leadDiffAxons.Saccadic],A(1,:),[],true,false);
subplot(1,2,2)
transform_swc_AV([lagDiffAxons.Saccadic],A(256,:),[],true,false);
export_fig('/Users/ashwin/Desktop/LeadLag_Saccadic_to_ABD_R_ABDi_B.png','-r300','-transparent');
close all;

figure('units','normalized','outerposition',[0 0 1 1]);
subplot(1,2,1)
transform_swc_AV(SaccadicOntoABD,A(128-64,:),[],true,false);
subplot(1,2,2)
transform_swc_AV(SaccadicOntoABDi,A(128+64,:),[],true,false);
clear temp;
export_fig('/Users/ashwin/Desktop/AllSaccadic_to_ABD_R_ABDi_B.png','-r300','-transparent');
close all;

% Vest pop
figure('units','normalized','outerposition',[0 0 1 1]);
temp = [sum(vestibularMotorNeuronTargets(:,2:3),2)-sum(vestibularMotorNeuronTargets(:,4:5),2)];
VestibularOntoABD = uniqueVestibularAxons(temp>0);
VestibularOntoABDi = uniqueVestibularAxons(temp<0);

subplot(1,2,1)
transform_swc_AV([leadDiffAxons.Vestibular],A(1,:),[],true,false);
subplot(1,2,2)
transform_swc_AV([lagDiffAxons.Vestibular],A(256,:),[],true,false);
export_fig('/Users/ashwin/Desktop/LeadLag_Vestibular_to_ABD_R_ABDi_B.png','-r300','-transparent');
close all;

figure('units','normalized','outerposition',[0 0 1 1]);
subplot(1,2,1)
transform_swc_AV(VestibularOntoABD,A(128-64,:),[],true,false);
subplot(1,2,2)
transform_swc_AV(VestibularOntoABDi,A(128+64,:),[],true,false);
clear temp;
export_fig('/Users/ashwin/Desktop/AllVestibular_to_ABD_R_ABDi_B.png','-r300','-transparent');
close all;

 % Integrator Pop
figure('units','normalized','outerposition',[0 0 1 1]);
temp = [sum(integratorMotorNeuronTargets(:,2:3),2)-sum(integratorMotorNeuronTargets(:,4:5),2)];
IntegratorOntoABD = uniqueIntegratorAxons(temp>0);
IntegratorOntoABDi = uniqueIntegratorAxons(temp<0);

subplot(1,2,1)
transform_swc_AV([leadDiffAxons.Integrator],A(1,:),[],true,false);
subplot(1,2,2)
transform_swc_AV([lagDiffAxons.Integrator],A(256,:),[],true,false);
export_fig('/Users/ashwin/Desktop/LeadLag_Integrator_to_ABD_R_ABDi_B.png','-r300','-transparent');
close all;

figure('units','normalized','outerposition',[0 0 1 1]);
subplot(1,2,1)
transform_swc_AV(IntegratorOntoABD,A(128-64,:),[],true,false);
subplot(1,2,2)
transform_swc_AV(IntegratorOntoABDi,A(128+64,:),[],true,false);
clear temp;
export_fig('/Users/ashwin/Desktop/AllIntegrator_to_ABD_R_ABDi_B.png','-r300','-transparent');
close all;

%%

% Contra pop
figure('units','normalized','outerposition',[0 0 1 1]);
temp = [sum(contraMotorNeuronTargets(:,2:3),2)-sum(contraMotorNeuronTargets(:,4:5),2)];
ContraOntoABD = uniqueContraAxons(temp>0);
ContraOntoABDi = uniqueContraAxons(temp<0);

subplot(1,2,1)
transform_swc_AV([leadDiffAxons.Contra],A(1,:),[],true,false);
subplot(1,2,2)
transform_swc_AV([lagDiffAxons.Contra],A(256,:),[],true,false);
export_fig('/Users/ashwin/Desktop/LeadLag_Contra_to_ABD_R_ABDi_B.png','-r300','-transparent');
close all;

figure('units','normalized','outerposition',[0 0 1 1]);
subplot(1,2,1)
transform_swc_AV(ContraOntoABD,A(128-64,:),[],true,false);
subplot(1,2,2)
transform_swc_AV(ContraOntoABDi,A(128+64,:),[],true,false);
clear temp;
export_fig('/Users/ashwin/Desktop/AllContra_to_ABD_R_ABDi_B.png','-r300','-transparent');
close all;



%%
figure('units','normalized','outerposition',[0 0 1 1]);
temp = [sum(EverythingElseAxonsMotorNeuronTargets(:,2:3),2)-sum(EverythingElseAxonsMotorNeuronTargets(:,4:5),2)];
EverythingElseOntoABD = uniqueEverythingElseAxons(temp>0);
EverythingElseOntoABDi = uniqueEverythingElseAxons(temp<0);
transform_swc_AV(IntegratorOntoABD,A(1,:),[],false);
transform_swc_AV(IntegratorOntoABDi,A(3,:),[],false);
clear temp


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







