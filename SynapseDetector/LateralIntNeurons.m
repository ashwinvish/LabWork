
% lateral class

clear;
addpath(genpath('/Users/ashwin/Documents/'));
colorSchemes;

startup

if ismac
    addpath(genpath('/Users/ashwin/Documents/'));
    df = readtable('/Users/ashwin/Google Drive/Zfish/SynapseDetector/04152019.csv');
else
    addpath(genpath('/usr/people/ashwinv/seungmount/research/Ashwin/Scripts'));
    df = readtable('/usr/people/ashwinv/seungmount/research/Ashwin/SynapseDetector/11252018.csv');
end



lateralVSaccadic.cellIDs = [80163 80167 80177 80179 76688 80204 80206 80210 76682];
 
vestibularCellIds = [76631,76700,76656,76655,76660,76679,77393,77395,77375,...
     77815,77803,77807,80591,80579,77255,77697,77364,80223,81126,80760,81117,...
     80672,80622,80572,80569,81122,78426,81328,81575,81582,81870,82161,82165,81553,82169];
 
 MVNs = [76664 76784 76783 77394 77326 77772 77766 77756 77753 77767 77775 ...
     77780 80883 80247 80669 80698 80762 80748 80696 80761 80749 81070 81044 ...
     80991 80967 81008 80883 81045 80986 81023 78647 78619 80247 81141];
 
VestibularAxons = [vestibularCellIds,MVNs];

load ABDr.mat
load ABDc.mat



% transform_swc_AV(lateralVSaccadic.cellIDs,SaccABDcolor,[],true,true,'LatInt-ABD-VEL');
% transform_swc_AV(lateralVSaccadic.cellIDs,SaccABDcolor,[],true,true,'LatInt-ABD-POS');

% output gradients

[lateralVSaccadic.ABDrGradientPathLength,lateralVSaccadic.ABDrGradient] = getABDgradient('ABDr',lateralVSaccadic.cellIDs',true,false);
[lateralVSaccadic.ABDcGradientPathLength,lateralVSaccadic.ABDcGradient] = getABDgradient('ABDc',lateralVSaccadic.cellIDs',true,false);

figure;
subplot(4,4,1)
errorbar(nanmean([lateralVSaccadic.ABDrGradient;lateralVSaccadic.ABDcGradient]),nanstd([lateralVSaccadic.ABDrGradient;lateralVSaccadic.ABDcGradient])./sqrt(29),...
    '-o','color','r','LineWidth',2,'MarkerFaceColor','w');
xlabel('Norm. pathlength');
ylabel('Norm. count');
axis square;
box off;
offsetAxes(gca);

% input gradients

for i = 1:length(lateralVSaccadic.cellIDs)
    Lat(i) = InputsByClass(lateralVSaccadic.cellIDs(i),df,1);
end


load ABDPutativeSaccadic.mat
load ABDiPutativeSaccadic.mat
load SaccadicProjectingToABDexclusively.mat

[lateralVSaccadic.r23ABDPathLenght,lateralVSaccadic.r23ABDGradient] = getABDgradient(Lat,ABDPutativeSaccadic.cellIDs',true,false);
[lateralVSaccadic.r23ABDiPathLenght,lateralVSaccadic.r23ABDiGradient] = getABDgradient(Lat,ABDiPutativeSaccadic.cellIDs',true,false);
[lateralVSaccadic.VestibularPathLength,lateralVSaccadic.VestibularGradient] = getABDgradient(Lat,VestibularAxons',true,false);
[lateralVSaccadic.r456ABDPathLenght,lateralVSaccadic.r456ABDGradient] = getABDgradient(Lat,SaccadicProjectingToABDexclusively.cellID',true,false);

figure;
subplot(1,4,1)
heatmap(lateralVSaccadic.r23ABDGradient);
title('Sacc');

subplot(1,4,2)
heatmap(lateralVSaccadic.r23ABDiGradient);
title('Sacci');

subplot(1,4,3)
heatmap(lateralVSaccadic.VestibularGradient);
title('Vest');

subplot(1,4,4)
heatmap(lateralVSaccadic.r456ABDGradient);
title('r456');

figure;

subplot(4,4,1)
histogram([cell2mat(lateralVSaccadic.r23ABDPathLenght)';cell2mat(lateralVSaccadic.r23ABDiPathLenght)'],20,'Normalization','cdf',...
    'DisplayStyle','stairs','LineWidth',2,'EdgeColor','k');
hold on;
histogram(lateralVSaccadic.VestibularGradient,20,'Normalization','cdf',...
    'DisplayStyle','stairs','LineWidth',2,'EdgeColor','r');
histogram(lateralVSaccadic.r456ABDGradient,20,'Normalization','cdf',...
    'DisplayStyle','stairs','LineWidth',2,'EdgeColor','b');
box off;
axis square;
set(gca, 'XLim',[0,1]);
xlabel('Norm.pathlength');
ylabel('Norm. count');
offsetAxes(gca);





%%
load ABDPutativeSaccadic.mat
load ABDiPutativeSaccadic.mat
load ABD_ABDi_PutativeSaccadic.mat

ABDPutativeSaccadic.postIBNs = [];

for i = 1:size(ABDPutativeSaccadic.cellIDs,1)
temp = SynapticPartners(ABDPutativeSaccadic.cellIDs(i),2,df);
temp = temp(temp<1e5);
ABDPutativeSaccadic.postIBNs = [ABDPutativeSaccadic.postIBNs;temp(isIBN(temp))];
clear temp;
end

ABDiPutativeSaccadic.postIBNs = [];

for i = 1:size(ABDiPutativeSaccadic.cellIDs,1)
temp = SynapticPartners(ABDiPutativeSaccadic.cellIDs(i),2,df);
temp = temp(temp<1e5);
ABDiPutativeSaccadic.postIBNs = [ABDiPutativeSaccadic.postIBNs;temp(isIBN(temp))];
clear temp;
end

%% plot vel Pos maps for IBNs
transform_swc_AV(IBNall,IbnABDcolor,[],true,true,'IBNall-Vel');
transform_swc_AV(IBNall,IbnABDcolor,[],true,true,'IBNall-POS');


