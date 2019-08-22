% Location of Saccadic synapses onto IBNs

startup

if ismac
    addpath(genpath('/Users/ashwin/Documents/'));
    df = readtable('/Users/ashwin/Google Drive/Zfish/SynapseDetector/04152019.csv');
else
    addpath(genpath('/usr/people/ashwinv/seungmount/research/Ashwin/Scripts'));
    df = readtable('/usr/people/ashwinv/seungmount/research/Ashwin/SynapseDetector/11252018.csv');
end

%IBNordered = [77125 77128 77231 77941 77940 77157 77247 77153 79053 77135 77941 78550];
%IBNrem = [77137 77942 78557 78567 78685 79083 79084 80971];

IBNall = [ 77125 77128 77135 77153 77231 77247 78685 77941 77137 79053 77940 77942 80971 77157 78550 79084 78557 79083 78567];

%IBNall = [SaccadicToIBN.ABDcellIDs,SaccadicToIBN.ABDicellIDs];
%IBNall = [IBNordered,IBNrem];

for i = 1:numel(IBNall)
    IBN(i) = InputsByClass(IBNall(i),df);
end


%%
% Plot Saccadic distributions

for i = 1:numel(IBNall)
    if ~isempty(IBN(i).Tree)
    IBN(i).SaccadicDist  = IBN(i).PathLength(IBN(i).isSaccadic)/max(Pvec_tree(IBN(i).Tree{1}));
    IBN(i).VestibularDist  = IBN(i).PathLength(IBN(i).isVestibular)/max(Pvec_tree(IBN(i).Tree{1}));
    IBN(i).ContraDist  = IBN(i).PathLength(IBN(i).isContra)/max(Pvec_tree(IBN(i).Tree{1}));
    IBN(i).IntegratorDist  = IBN(i).PathLength(IBN(i).isIntegrator)/max(Pvec_tree(IBN(i).Tree{1}));
    IBN(i).EverythingElseDist = IBN(i).PathLength(IBN(i).isEverythingElse)/max(Pvec_tree(IBN(i).Tree{1}));
    end
end

%% plot distributions based on ordering

subplot(4,4,1)
histogram(vertcat(IBN(1:6).SaccadicDist),20,'DisplayStyle','stairs','EdgeColor',IbnABDcolor,'LineWidth',2);
hold on
histogram(vertcat(IBN(7:end).SaccadicDist),20,'DisplayStyle','stairs','EdgeColor',IbnABDicolor,'LineWidth',2);
axis square;
box off;
legend('Sm->Im','Si-->Ii','Location','best');
xlabel('Norm. pathlength');
ylabel('Count');
offsetAxes(gca);


subplot(4,4,2)
histogram(vertcat(IBN(1:6).VestibularDist),20,'DisplayStyle','stairs','EdgeColor',IbnABDcolor,'LineWidth',2);
hold on
histogram(vertcat(IBN(7:end).VestibularDist),20,'DisplayStyle','stairs','EdgeColor',IbnABDicolor,'LineWidth',2)
axis square;
box off;
%legend('Sm->Im','Si-->Ii');

subplot(4,4,3)
histogram(vertcat(IBN(1:6).ContraDist),20,'DisplayStyle','stairs','EdgeColor',IbnABDcolor,'LineWidth',2);
hold on
histogram(vertcat(IBN(7:end).ContraDist),20,'DisplayStyle','stairs','EdgeColor',IbnABDicolor,'LineWidth',2)
axis square;
box off;


%%
subplot(4,5,1)
h1 = histogram(vertcat(IBN.SaccadicDist),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = [1,0,0];
h1.LineWidth = 2;
set(gca,'XLim',[0,1],'YLim',[0,0.15]);
box off;

subplot(4,5,2)
h1 = histogram(vertcat(IBN.VestibularDist),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = [1,0,0];
h1.LineWidth = 2;
set(gca,'XLim',[0,1],'YLim',[0,0.15]);
box off;

subplot(4,5,3)
h1 = histogram(vertcat(IBN.IntegratorDist),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = [1,0,0];
h1.LineWidth = 2;
set(gca,'XLim',[0,1],'YLim',[0,0.15]);
box off;

subplot(4,5,4)
h1 = histogram(vertcat(IBN.ContraDist),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = [1,0,0];
h1.LineWidth = 2;
set(gca,'XLim',[0,1],'YLim',[0,0.15]);
box off;

subplot(4,5,5)
h1 = histogram(vertcat(IBN.EverythingElseDist),20,'Normalization','probability','FaceColor','none','DisplayStyle','stairs');
h1.EdgeColor = [1,0,0];
h1.LineWidth = 2;
set(gca,'XLim',[0,1],'YLim',[0,0.15]);
box off;


%% check for IBN contra pop.