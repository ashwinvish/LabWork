% Location of Saccadic synapses onto IBNs



if ismac
    addpath(genpath('/Users/ashwin/Documents/'));
    df = readtable('/Users/ashwin/Google Drive/Zfish/SynapseDetector/04152019.csv');
else
    addpath(genpath('/usr/people/ashwinv/seungmount/research/Ashwin/Scripts'));
    df = readtable('/usr/people/ashwinv/seungmount/research/Ashwin/SynapseDetector/11252018.csv');
end
startup

load ABDPutativeSaccadic.mat
load ABDiPutativeSaccadic.mat

lateralVSaccadic = [80163 80167 80177 80179 76688 80204 80206 80210 76682];

%IBNordered = [77125 77128 77231 77941 77940 77157 77247 77153 79053 77135 77941 78550];
%IBNrem = [77137 77942 78557 78567 78685 79083 79084 80971];

IBNall = [77153,77231,82710,77137,82709,77125,82713,78685,77942,82264,77247,83184,78550,83183,82711];

%IBNall = [SaccadicToIBN.Sac_ABDcellIDs, SaccadicToIBN.Sac_ABDicellIDs];
%IBNall = [IBNordered,IBNrem];

for i = 1:numel(IBNall)
    IBN(i) = InputsByClass(IBNall(i),df,1);
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
    temp = find(ismember(IBN(i).Inputs,ABDPutativeSaccadic.cellIDs'));
    IBN(i).ABDputativeSacc = IBN(i).PathLength(temp)/max(Pvec_tree(IBN(i).Tree{1}));
    temp2 = find(ismember(IBN(i).Inputs,ABDiPutativeSaccadic.cellIDs'));
    IBN(i).ABDiputativeSacc = IBN(i).PathLength(temp2)/max(Pvec_tree(IBN(i).Tree{1}));
    temp3 = find(ismember(IBN(i).Inputs,lateralVSaccadic));
     IBN(i).r5integ = IBN(i).PathLength(temp3)/max(Pvec_tree(IBN(i).Tree{1}));
    end
    clear temp
    clear temp2
    clear temp3
end

subplot(4,4,1)
 histogram(vertcat(IBN.ABDputativeSacc,IBN.ABDiputativeSacc),20,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'Edgecolor',[74,74,74]./255);
 hold on;
 histogram(vertcat(IBN.VestibularDist),20,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'Edgecolor',[255,127,0]./255);
 histogram(vertcat(IBN.IntegratorDist,IBN.SaccadicDist,IBN.r5integ),20,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'Edgecolor',[55,126,184]./255);

axis square;
box off;
offsetAxes(gca);
xlabel('Norm. pathlength');
ylabel('cdf');
legend({'Sacc.','Vest.','Integ.'},'Location','bestoutside');

%% gradients
saccadicColorMap = ['#ffffff'; '#f5f5f5'; '#eaeaea'; '#e0e0e0'; '#d6d6d6'; '#cccccc'; '#c2c2c2'; '#b8b8b8'; '#aeaeae'; '#a4a4a4'; '#9b9b9b'; '#919191'; '#888888'; '#7e7e7e'; '#757575'; '#6c6c6c'; '#636363'; '#5b5b5b'; '#525252'; '#4a4a4a'];
%saccadicColorMap = colorcet('L1','N',20,'reverse',1);
vestibuarColorMap = ['#ffffff'; '#fff9f4'; '#fff3e9'; '#ffedde'; '#ffe7d3'; '#ffe1c7'; '#ffdabc'; '#ffd4b1'; '#ffcea6'; '#ffc79a'; '#ffc18e'; '#ffba83'; '#ffb477'; '#ffad6a'; '#ffa65e'; '#ff9e50'; '#ff9742'; '#ff8f33'; '#ff8720'; '#ff7f00'];
integratorColorMap = ['#ffffff'; '#f5f8fb'; '#ecf1f8'; '#e2e9f4'; '#d9e2f0'; '#cfdbec'; '#c5d4e9'; '#bccde5'; '#b2c6e1'; '#a8c0dd'; '#9eb9da'; '#94b2d6'; '#8aabd2'; '#80a5ce'; '#769ecb'; '#6b98c7'; '#5f91c3'; '#538bbf'; '#4684bc'; '#377eb8'];

subplot(2,3,1)
[a,b] = getABDgradient(IBN,[ABDPutativeSaccadic.cellIDs';ABDiPutativeSaccadic.cellIDs'],false);
heatmap(b,'Colormap',hex2rgb(saccadicColorMap),'ColorbarVisible','on','XDisplayLabels',0.1:0.1:1);
title('r2/3 Saccadic');


subplot(2,3,2)

[c,d] = getABDgradient(IBN,unique(vertcat(IBN.Vestibular)),false);
heatmap(d,'Colormap',hex2rgb(vestibuarColorMap),'ColorbarVisible','on','XDisplayLabels',0.1:0.1:1);
title('vestibular');


subplot(2,3,3)

[e,f] = getABDgradient(IBN,[unique(vertcat(IBN.Saccadic));unique(vertcat(IBN.Integrator));lateralVSaccadic'],false);
heatmap(f,'Colormap',hex2rgb(integratorColorMap),'ColorbarVisible','on','XDisplayLabels',0.1:0.1:1);
title('Integrator');

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

% fraction reconstructed

for i = 1:numel(IBNall)
    IBN(i).ReconFraction = sum(IBN(i).Inputs<1e5)/numel(IBN(i).Inputs);
end

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

%% IBN manual contra pop.
 
manualPop = [77303,79799,80194,80847,80701,77355,77353,80737];
manualPop_motorDist = isMotor(manualPop',df);

%manualPop = [80816,80707,81169,77855,79799,80847,77355,80700,80868];

%% check for IBN contra pop.

IBNSacc.ContraCellIDs = [];
for i = 1:numel(SaccadicToIBN.Sac_ABDcellIDs)
    IBNSacc.ContraCellIDs = [IBNSacc.ContraCellIDs;IBN(i).Inputs(IBN(i).isContra)];
end
 IBNSacc.ContraCellIDs = unique( IBNSacc.ContraCellIDs);
 IBNSacc.ContraRhombomere = isRhombomere(IBNSacc.ContraCellIDs);

IBNSacci.ContraCellIDs = [];
for i = 1:numel(SaccadicToIBN.Sac_ABDicellIDs)
    IBNSacci.ContraCellIDs = [IBNSacci.ContraCellIDs;IBN(i).Inputs(IBN(i).isContra)];
end

IBNSacci.ContraCellIDs = unique( IBNSacci.ContraCellIDs);
IBNSacci.ContraRhombomere = isRhombomere(IBNSacci.ContraCellIDs);

figure;
subplot(1,2,1)
 transform_swc_AV(IBNSacc.ContraCellIDs,IbnABDcolor,[],true,false);
 subplot(1,2,2)
 transform_swc_AV(IBNSacci.ContraCellIDs,IbnABDicolor,[],true,false);

 
 %% IBN gradients
 
[~,IBN.SaccadicGradient] = getABDgradient(IBN,SaccadicToIBN.Sac_ABDcellIDs,true);

