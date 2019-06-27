clear;
addpath(genpath('/Users/ashwin/Documents/'));

colorPallete = cbrewer('div','BrBG',5);

smallColor = colorPallete(1,:);
largeColor = colorPallete(5,:);

small = 'MIF';
large = 'SIF';

colorPallete = colorcet('D2','N',5);

lightGreen = colorPallete(1,:);
lightMagenta = colorPallete(5,:);


startup

if ismac
    addpath(genpath('/Users/ashwin/Documents/'));
    df = readtable('/Users/ashwin/Documents/SynapseDetector/04152019.csv');
else
    addpath(genpath('/usr/people/ashwinv/seungmount/research/Ashwin/Scripts'));
    df = readtable('/usr/people/ashwinv/seungmount/research/Ashwin/SynapseDetector/11252018.csv');
end

ABDr_CellIDs = [82140,82145,82143,77648,82146,77302,77301, 82194,77705,77710,77305,77672,77300,77709];
ABDc_CellIDs = [82213,77654,77646,82212,77295,81172,77657,82197,77682,77154, 77652,77658,77628,77292,77688,82195,77296];
ABDIr_CellIDs = [78553,77150,77886,78547,77631,77158,77665,78552,77668,77618,77634];
ABDIc_CellIDs = [78574,79051,77148,77625,77144,77692,77641,79066,77640,77643];


AllABD = [ABDr_CellIDs,ABDc_CellIDs,ABDIr_CellIDs,ABDIc_CellIDs];
Allmotor = [ABDr_CellIDs,ABDc_CellIDs];
%%
load('ABDVols.mat');

for i = 1:size(ABDvols,1)
    ABD.vol(i,1) = ABDvols(i,2);
    ABD.cellID(i,1) = ABDvols(i,1);
    ABD.Synapses(i,1) = size(SynapticPartners(ABDvols(i,1),1,df),1);
    if (isExistReRoot(ABDvols(i,1)) == 1)
        tempTree = SwctoZbrian(ABDvols(i,1));
        ABD.PathLength(i,1) = max(Pvec_tree(tempTree{1}));
        ABD.Origin(i,:) = getOrigin(ABDvols(i,1));
        clear tempTree;
    else
        ABD.PathLength(i,1)= NaN;
        ABD.Origin(i,:) = [NaN,NaN NaN];
    end
end

[~,ABD.sortOrder] = sort(ABDvols(:,2));
ABD.vol = ABD.vol/1e9;


[~,~,ABD.ABDrLocs] = intersect(ABDr_CellIDs,ABD.cellID,'stable');
ABD.ABDrvols = ABD.vol(ABD.ABDrLocs);
ABD.ABDrTotalSynapses = ABD.Synapses(ABD.ABDrLocs);

[~,~,ABD.ABDcLocs] = intersect(ABDc_CellIDs,ABD.cellID,'stable');
ABD.ABDcvols = ABD.vol(ABD.ABDcLocs);
ABD.ABDcTotalSynapses = ABD.Synapses(ABD.ABDcLocs);

[~,~,ABD.ABDIrLocs] = intersect(ABDIr_CellIDs,ABD.cellID,'stable');
ABD.ABDIrvols = ABD.vol(ABD.ABDIrLocs);
ABD.ABDIrTotalSynapses =  ABD.Synapses(ABD.ABDIrLocs);

[~,~,ABD.ABDIcLocs] = intersect(ABDIc_CellIDs,ABD.cellID,'stable');
ABD.ABDIcvols = ABD.vol(ABD.ABDIcLocs);
ABD.ABDIcTotalSynapses =  ABD.Synapses(ABD.ABDIcLocs);

figure;
subplot(4,4,1);
scatter(ABD.Origin(ABD.ABDcLocs), ABD.ABDcvols,25,lightGreen,'filled')
hold on;
scatter(ABD.Origin(ABD.ABDrLocs), ABD.ABDrvols,25,lightGreen,'filled')
scatter(ABD.Origin(ABD.ABDIrLocs), ABD.ABDIrvols,25,lightMagenta,'filled')
scatter(ABD.Origin(ABD.ABDIcLocs), ABD.ABDIcvols,25,lightMagenta,'filled')
axis square
xlabel('Medial <---> Lateral');
ylabel('Soma volume \mu^3');

subplot(4,4,2);
scatter(ABD.ABDrvols, ABD.ABDrTotalSynapses,25,lightGreen,'filled')
hold on;
scatter(ABD.ABDcvols, ABD.ABDcTotalSynapses,25,lightGreen,'filled')
scatter(ABD.ABDIrvols, ABD.ABDIrTotalSynapses,25,lightMagenta,'filled')
scatter(ABD.ABDIcvols,ABD.ABDIcTotalSynapses,25,lightMagenta,'filled')
axis square;
ylabel('Synapses');
xlabel('Soma volume \mu^3');

subplot(4,4,4)
scatter3([ABD.Origin(ABD.ABDrLocs);ABD.Origin(ABD.ABDcLocs)], [ABD.Origin(ABD.ABDrLocs,3);ABD.Origin(ABD.ABDcLocs,3)], [ABD.ABDrvols;ABD.ABDcvols],25,lightGreen,'filled')
hold on
scatter3([ABD.Origin(ABD.ABDIrLocs);ABD.Origin(ABD.ABDIcLocs)], [ABD.Origin(ABD.ABDIrLocs,3);ABD.Origin(ABD.ABDIcLocs,3)], [ABD.ABDIrvols;ABD.ABDIcvols],25,lightMagenta,'filled')
axis square;
ylabel('Dorsal <---> Ventral');
xlabel('Medial <---> Lateral');
zlabel('Soma volume \mu^3');
set(gca, 'YDir','reverse');

%%
load ABDr.mat
load ABDc.mat
load ABDIr.mat
load ABDIc.mat
%%

for i = 1:numel(ABDr)
    ABD.SaccadicSynapses(i) = length(ABDr(i).Saccadic);
    ABD.ABDrVestibularSynapses(i) = length(ABDr(i).Vestibular);
    ABD.ABDrIntegratorSynapses(i) = length(ABDr(i).Integrator);
    ABD.ABDrContraSynapses(i) = length(ABDr(i).Contra);
    ABD.ABDrRestSynapses(i) = length(ABDr(i).EverythingElse);
end

for i = 1:numel(ABDc)
    ABD.ABDcSaccadicSynapses(i) = length(ABDc(i).Saccadic);
    ABD.ABDcVestibularSynapses(i) = length(ABDc(i).Vestibular);
    ABD.ABDcIntegratorSynapses(i) = length(ABDc(i).Integrator);
    ABD.ABDcContraSynapses(i) = length(ABDc(i).Contra);
    ABD.ABDcRestSynapses(i) = length(ABDc(i).EverythingElse);
end


cols = cbrewer('qual','Dark2',5);
figure;
subplot(4,4,1);
scatter([ABD.ABDrvols;ABD.ABDcvols],[ABD.ABDrSaccadicSynapses./ABD.ABDrTotalSynapses',ABD.ABDcSaccadicSynapses./ABD.ABDcTotalSynapses'],25,cols(1,:),'filled','MarkerEdgeColor','k');
hold on;
showfit(ezfit([ABD.ABDrvols;ABD.ABDcvols],[ABD.ABDrSaccadicSynapses./ABD.ABDrTotalSynapses',ABD.ABDcSaccadicSynapses./ABD.ABDcTotalSynapses'],'lin'),'fitcolor',cols(1,:),'dispeqboxmode','off');
axis square;
xlabel('Soma volume \mu^3');
ylabel('Fraction of total synapses');


showfit(ezfit(ABD.ABDrvols,ABD.ABDrVestibularSynapses./ABD.ABDrTotalSynapses','affine'),'fitcolor',cols(2,:),'dispeqboxmode','off')
showfit(ezfit(ABD.ABDrvols,ABD.ABDrIntegratorSynapses./ABD.ABDrTotalSynapses','affine'),'fitcolor',cols(3,:),'dispeqboxmode','off')
showfit(ezfit(ABD.ABDrvols,ABD.ABDrContraSynapses./ABD.ABDrTotalSynapses','affine'),'fitcolor',cols(4,:),'dispeqboxmode','off')
showfit(ezfit(ABD.ABDrvols,ABD.ABDrRestSynapses./ABD.ABDrTotalSynapses','affine'),'fitcolor',cols(5,:),'dispeqboxmode','off')

    


%%
for i = 1:numel(ABDr_CellIDs)
    temp1 = ismember(ABDr(i).Inputs,unique([vertcat(ABDr.Saccadic);vertcat(ABDc.Saccadic);vertcat(ABDIr.Saccadic);vertcat(ABDIc.Saccadic)]));
    ABD.SaccPathLengthABDr{i} = ABDr(i).PathLength(temp1)'./max(Pvec_tree(ABDr(i).Tree{1}));
    clear temp1;
    
    temp1 = ismember(ABDr(i).Inputs,unique([vertcat(ABDr.Vestibular);vertcat(ABDc.Vestibular);vertcat(ABDIr.Vestibular);vertcat(ABDIc.Vestibular)]));
    ABD.VestPathLengthABDr{i} = ABDr(i).PathLength(temp1)'./max(Pvec_tree(ABDr(i).Tree{1}));
    clear temp1;
    
    temp1 = ismember(ABDr(i).Inputs,unique([vertcat(ABDr.Integrator);vertcat(ABDc.Integrator);vertcat(ABDIr.Integrator);vertcat(ABDIc.Integrator)]));
    ABD.IntPathLengthABDr{i} = ABDr(i).PathLength(temp1)'./max(Pvec_tree(ABDr(i).Tree{1}));
    clear temp1;
    
    temp1 = ismember(ABDr(i).Inputs,unique([vertcat(ABDr.Contra);vertcat(ABDc.Contra);vertcat(ABDIr.Contra);vertcat(ABDIc.Contra)]));
    ABD.ContraPathLengthABDr{i} = ABDr(i).PathLength(temp1)'./max(Pvec_tree(ABDr(i).Tree{1}));
    clear temp1
    
    temp1 = ismember(ABDr(i).Inputs,unique([vertcat(ABDr.EverythingElse);vertcat(ABDc.EverythingElse);vertcat(ABDIr.EverythingElse);vertcat(ABDIc.EverythingElse)]));
    ABD.RestPathLengthABDr{i} = ABDr(i).PathLength(temp1)'./max(Pvec_tree(ABDr(i).Tree{1}));
    clear temp1
end


for i = 1:numel(ABDr_CellIDs)
    ABD.ABDrSaccadicGradient(i,:) =  histcounts(ABD.SaccPathLengthABDr{i},0:0.1:1,'Normalization','probability')/max(histcounts(ABD.SaccPathLengthABDr{i},0:0.1:1,'Normalization','probability'));
    ABD.ABDrVestibularGradient(i,:) =  histcounts(ABD.VestPathLengthABDr{i},0:0.1:1,'Normalization','probability')/max(histcounts(ABD.VestPathLengthABDr{i},0:0.1:1,'Normalization','probability'));
    ABD.ABDrIntegratorGradient(i,:) =  histcounts(ABD.IntPathLengthABDr{i},0:0.1:1,'Normalization','probability')/max(histcounts(ABD.IntPathLengthABDr{i},0:0.1:1,'Normalization','probability'));
    ABD.ABDrContraGradient(i,:) =  histcounts(ABD.ContraPathLengthABDr{i},0:0.1:1,'Normalization','probability')/max(histcounts(ABD.ContraPathLengthABDr{i},0:0.1:1,'Normalization','probability'));
    ABD.ABDrRestGradient(i,:) =  histcounts(ABD.RestPathLengthABDr{i},0:0.1:1,'Normalization','probability')/max( histcounts(ABD.RestPathLengthABDr{i},0:0.1:1,'Normalization','probability'));
end

% ABDc

for i = 1:numel(ABDc_CellIDs)
    if ~isempty(ABDc(i).Tree)
        
        temp1 = ismember(ABDc(i).Inputs,unique([vertcat(ABDc.Saccadic);vertcat(ABDc.Saccadic);vertcat(ABDIr.Saccadic);vertcat(ABDIc.Saccadic)]));
        ABD.SaccPathLengthABDc{i} = ABDc(i).PathLength(temp1)'./max(Pvec_tree(ABDc(i).Tree{1}));
        clear temp1;
        
        temp1 = ismember(ABDc(i).Inputs,unique([vertcat(ABDc.Vestibular);vertcat(ABDc.Vestibular);vertcat(ABDIr.Vestibular);vertcat(ABDIc.Vestibular)]));
        ABD.VestPathLengthABDc{i} = ABDc(i).PathLength(temp1)'./max(Pvec_tree(ABDc(i).Tree{1}));
        clear temp1;
        
        temp1 = ismember(ABDc(i).Inputs,unique([vertcat(ABDc.Integrator);vertcat(ABDc.Integrator);vertcat(ABDIr.Integrator);vertcat(ABDIc.Integrator)]));
        ABD.IntPathLengthABDc{i} = ABDc(i).PathLength(temp1)'./max(Pvec_tree(ABDc(i).Tree{1}));
        clear temp1;
        
        temp1 = ismember(ABDc(i).Inputs,unique([vertcat(ABDc.Contra);vertcat(ABDc.Contra);vertcat(ABDIr.Contra);vertcat(ABDIc.Contra)]));
        ABD.ContraPathLengthABDc{i} = ABDc(i).PathLength(temp1)'./max(Pvec_tree(ABDc(i).Tree{1}));
        clear temp1
        
        temp1 = ismember(ABDc(i).Inputs,unique([vertcat(ABDc.EverythingElse);vertcat(ABDc.EverythingElse);vertcat(ABDIr.EverythingElse);vertcat(ABDIc.EverythingElse)]));
        ABD.RestPathLengthABDc{i} = ABDc(i).PathLength(temp1)'./max(Pvec_tree(ABDc(i).Tree{1}));
        clear temp1
    end
end


for i = 1:numel(ABDc_CellIDs)
    ABD.ABDcSaccadicGradient(i,:) =  histcounts(ABD.SaccPathLengthABDc{i},0:0.1:1,'Normalization','probability')/max(histcounts(ABD.SaccPathLengthABDc{i},0:0.1:1,'Normalization','probability'));
    ABD.ABDcVestibularGradient(i,:) =  histcounts(ABD.VestPathLengthABDc{i},0:0.1:1,'Normalization','probability')/max(histcounts(ABD.VestPathLengthABDc{i},0:0.1:1,'Normalization','probability'));
    ABD.ABDcIntegratorGradient(i,:) =  histcounts(ABD.IntPathLengthABDc{i},0:0.1:1,'Normalization','probability')/max(histcounts(ABD.IntPathLengthABDc{i},0:0.1:1,'Normalization','probability'));
    ABD.ABDcContraGradient(i,:) =  histcounts(ABD.ContraPathLengthABDc{i},0:0.1:1,'Normalization','probability')/max(histcounts(ABD.ContraPathLengthABDc{i},0:0.1:1,'Normalization','probability'));
    ABD.ABDcRestGradient(i,:) =  histcounts(ABD.RestPathLengthABDc{i},0:0.1:1,'Normalization','probability')/max(histcounts(ABD.RestPathLengthABDc{i},0:0.1:1,'Normalization','probability'));
end

% combine ABDr and ABDc
ABD.motorCellIDsSorted = intersect(ABD.sortOrder,find(ismember(ABD.cellID,Allmotor)),'stable');
ABD.motorVols = ABD.vol(ABD.motorCellIDsSorted);
ABD.motorSaccadicGradient(find(ismember(ABDvols(:,1),ABDr_CellIDs)),:) = ABD.ABDrSaccadicGradient;
ABD.motorSaccadicGradient(find(ismember(ABDvols(:,1),ABDc_CellIDs)),:) = ABD.ABDcSaccadicGradient;
ABD.motorSaccadicGradient(setdiff(1:52,[find(ismember(ABDvols(:,1),ABDr_CellIDs));find(ismember(ABDvols(:,1),ABDc_CellIDs))])',:) = [];

ABD.motorVestibularGradient(find(ismember(ABDvols(:,1),ABDr_CellIDs)),:) = ABD.ABDrVestibularGradient;
ABD.motorVestibularGradient(find(ismember(ABDvols(:,1),ABDc_CellIDs)),:) = ABD.ABDcVestibularGradient;
ABD.motorVestibularGradient(setdiff(1:42,[find(ismember(ABDvols(:,1),ABDr_CellIDs));find(ismember(ABDvols(:,1),ABDc_CellIDs))])',:) = [];

ABD.motorIntegratorGradient(find(ismember(ABDvols(:,1),ABDr_CellIDs)),:) = ABD.ABDrIntegratorGradient;
ABD.motorIntegratorGradient(find(ismember(ABDvols(:,1),ABDc_CellIDs)),:) = ABD.ABDcIntegratorGradient;
ABD.motorIntegratorGradient(setdiff(1:42,[find(ismember(ABDvols(:,1),ABDr_CellIDs));find(ismember(ABDvols(:,1),ABDc_CellIDs))])',:) = [];

ABD.motorContraGradient(find(ismember(ABDvols(:,1),ABDr_CellIDs)),:) = ABD.ABDrContraGradient;
ABD.motorContraGradient(find(ismember(ABDvols(:,1),ABDc_CellIDs)),:) = ABD.ABDcContraGradient;
ABD.motorContraGradient(setdiff(1:42,[find(ismember(ABDvols(:,1),ABDr_CellIDs));find(ismember(ABDvols(:,1),ABDc_CellIDs))])',:) = [];

ABD.motorRestGradient(find(ismember(ABDvols(:,1),ABDr_CellIDs)),:) = ABD.ABDrRestGradient;
ABD.motorRestGradient(find(ismember(ABDvols(:,1),ABDc_CellIDs)),:) = ABD.ABDcRestGradient;
ABD.motorRestGradient(setdiff(1:42,[find(ismember(ABDvols(:,1),ABDr_CellIDs));find(ismember(ABDvols(:,1),ABDc_CellIDs))])',:) = [];


%%
% ABDIr

for i = 1:numel(ABDIr_CellIDs)
    temp1 = ismember(ABDIr(i).Inputs,unique([vertcat(ABDIr.Saccadic);vertcat(ABDIr.Saccadic);vertcat(ABDIr.Saccadic);vertcat(ABDIc.Saccadic)]));
    ABD.SaccPathLengthABDIr{i} = ABDIr(i).PathLength(temp1)'./max(Pvec_tree(ABDIr(i).Tree{1}));
    clear temp1;
    
    temp1 = ismember(ABDIr(i).Inputs,unique([vertcat(ABDIr.Vestibular);vertcat(ABDIr.Vestibular);vertcat(ABDIr.Vestibular);vertcat(ABDIc.Vestibular)]));
    ABD.VestPathLengthABDIr{i} = ABDIr(i).PathLength(temp1)'./max(Pvec_tree(ABDIr(i).Tree{1}));
    clear temp1;
    
    temp1 = ismember(ABDIr(i).Inputs,unique([vertcat(ABDIr.Integrator);vertcat(ABDIr.Integrator);vertcat(ABDIr.Integrator);vertcat(ABDIc.Integrator)]));
    ABD.IntPathLengthABDIr{i} = ABDIr(i).PathLength(temp1)'./max(Pvec_tree(ABDIr(i).Tree{1}));
    clear temp1;
    
    temp1 = ismember(ABDIr(i).Inputs,unique([vertcat(ABDIr.Contra);vertcat(ABDIr.Contra);vertcat(ABDIr.Contra);vertcat(ABDIc.Contra)]));
    ABD.ContraPathLengthABDIr{i} = ABDIr(i).PathLength(temp1)'./max(Pvec_tree(ABDIr(i).Tree{1}));
    clear temp1
    
    temp1 = ismember(ABDIr(i).Inputs,unique([vertcat(ABDIr.EverythingElse);vertcat(ABDIr.EverythingElse);vertcat(ABDIr.EverythingElse);vertcat(ABDIc.EverythingElse)]));
    ABD.RestPathLengthABDIr{i} = ABDIr(i).PathLength(temp1)'./max(Pvec_tree(ABDIr(i).Tree{1}));
    clear temp1
end


for i = 1:numel(ABDIr_CellIDs)
    ABD.ABDIrSaccadicGradient(i,:) =  histcounts(ABD.SaccPathLengthABDIr{i},0:0.1:1,'Normalization','probability')/max(histcounts(ABD.SaccPathLengthABDIr{i},0:0.1:1,'Normalization','probability'));
    ABD.ABDIrVestibularGradient(i,:) =  histcounts(ABD.VestPathLengthABDIr{i},0:0.1:1,'Normalization','probability')/max(histcounts(ABD.VestPathLengthABDIr{i},0:0.1:1,'Normalization','probability'));
    ABD.ABDIrIntegratorGradient(i,:) =  histcounts(ABD.IntPathLengthABDIr{i},0:0.1:1,'Normalization','probability')/max(histcounts(ABD.IntPathLengthABDIr{i},0:0.1:1,'Normalization','probability'));
    ABD.ABDIrContraGradient(i,:) =  histcounts(ABD.ContraPathLengthABDIr{i},0:0.1:1,'Normalization','probability')/max( histcounts(ABD.ContraPathLengthABDIr{i},0:0.1:1,'Normalization','probability'));
    ABD.ABDIrRestGradient(i,:) =  histcounts(ABD.RestPathLengthABDIr{i},0:0.1:1,'Normalization','probability')/max(histcounts(ABD.RestPathLengthABDIr{i},0:0.1:1,'Normalization','probability'));
end

% ABDIc

for i = 1:numel(ABDIc_CellIDs)
    temp1 = ismember(ABDIc(i).Inputs,unique([vertcat(ABDIc.Saccadic);vertcat(ABDIc.Saccadic);vertcat(ABDIc.Saccadic);vertcat(ABDIc.Saccadic)]));
    ABD.SaccPathLengthABDIc{i} = ABDIc(i).PathLength(temp1)'./max(Pvec_tree(ABDIc(i).Tree{1}));
    clear temp1;
    
    temp1 = ismember(ABDIc(i).Inputs,unique([vertcat(ABDIc.Vestibular);vertcat(ABDIc.Vestibular);vertcat(ABDIc.Vestibular);vertcat(ABDIc.Vestibular)]));
    ABD.VestPathLengthABDIc{i} = ABDIc(i).PathLength(temp1)'./max(Pvec_tree(ABDIc(i).Tree{1}));
    clear temp1;
    
    temp1 = ismember(ABDIc(i).Inputs,unique([vertcat(ABDIc.Integrator);vertcat(ABDIc.Integrator);vertcat(ABDIc.Integrator);vertcat(ABDIc.Integrator)]));
    ABD.IntPathLengthABDIc{i} = ABDIc(i).PathLength(temp1)'./max(Pvec_tree(ABDIc(i).Tree{1}));
    clear temp1;
    
    temp1 = ismember(ABDIc(i).Inputs,unique([vertcat(ABDIc.Contra);vertcat(ABDIc.Contra);vertcat(ABDIc.Contra);vertcat(ABDIc.Contra)]));
    ABD.ContraPathLengthABDIc{i} = ABDIc(i).PathLength(temp1)'./max(Pvec_tree(ABDIc(i).Tree{1}));
    clear temp1
    
    temp1 = ismember(ABDIc(i).Inputs,unique([vertcat(ABDIc.EverythingElse);vertcat(ABDIc.EverythingElse);vertcat(ABDIc.EverythingElse);vertcat(ABDIc.EverythingElse)]));
    ABD.RestPathLengthABDIc{i} = ABDIc(i).PathLength(temp1)'./max(Pvec_tree(ABDIc(i).Tree{1}));
    clear temp1
end


for i = 1:numel(ABDIc_CellIDs)
    ABD.ABDIcSaccadicGradient(i,:) =  histcounts(ABD.SaccPathLengthABDIc{i},0:0.1:1,'Normalization','probability')/max(histcounts(ABD.SaccPathLengthABDIc{i},0:0.1:1,'Normalization','probability'));
    ABD.ABDIcVestibularGradient(i,:) =  histcounts(ABD.VestPathLengthABDIc{i},0:0.1:1,'Normalization','probability')/max(histcounts(ABD.VestPathLengthABDIc{i},0:0.1:1,'Normalization','probability'));
    ABD.ABDIcIntegratorGradient(i,:) =  histcounts(ABD.IntPathLengthABDIc{i},0:0.1:1,'Normalization','probability')/max(histcounts(ABD.IntPathLengthABDIc{i},0:0.1:1,'Normalization','probability'));
    ABD.ABDIcContraGradient(i,:) =  histcounts(ABD.ContraPathLengthABDIc{i},0:0.1:1,'Normalization','probability')/max(histcounts(ABD.ContraPathLengthABDIc{i},0:0.1:1,'Normalization','probability'));
    ABD.ABDIcRestGradient(i,:) =  histcounts(ABD.RestPathLengthABDIc{i},0:0.1:1,'Normalization','probability')/max( histcounts(ABD.RestPathLengthABDIc{i},0:0.1:1,'Normalization','probability'));
end

% combine ABDIr and ABDIc
ABD.InterCellIDsSorted = intersect(ABD.sortOrder,find(setdiff(ABD.cellID,Allmotor)),'stable');
ABD.InterVols = ABD.vol(ABD.InterCellIDsSorted);
ABD.InterSaccadicGradient(find(ismember(ABDvols(:,1),ABDIr_CellIDs)),:) = ABD.ABDIrSaccadicGradient;
ABD.InterSaccadicGradient(find(ismember(ABDvols(:,1),ABDIc_CellIDs)),:) = ABD.ABDIcSaccadicGradient;
ABD.InterSaccadicGradient(setdiff(1:42,[find(ismember(ABDvols(:,1),ABDIr_CellIDs));find(ismember(ABDvols(:,1),ABDIc_CellIDs))])',:) = [];

ABD.InterVestibularGradient(find(ismember(ABDvols(:,1),ABDIr_CellIDs)),:) = ABD.ABDIrVestibularGradient;
ABD.InterVestibularGradient(find(ismember(ABDvols(:,1),ABDIc_CellIDs)),:) = ABD.ABDIcVestibularGradient;
ABD.InterVestibularGradient(setdiff(1:42,[find(ismember(ABDvols(:,1),ABDIr_CellIDs));find(ismember(ABDvols(:,1),ABDIc_CellIDs))])',:) = [];

ABD.InterIntegratorGradient(find(ismember(ABDvols(:,1),ABDIr_CellIDs)),:) = ABD.ABDIrIntegratorGradient;
ABD.InterIntegratorGradient(find(ismember(ABDvols(:,1),ABDIc_CellIDs)),:) = ABD.ABDIcIntegratorGradient;
ABD.InterIntegratorGradient(setdiff(1:42,[find(ismember(ABDvols(:,1),ABDIr_CellIDs));find(ismember(ABDvols(:,1),ABDIc_CellIDs))])',:) = [];

ABD.InterContraGradient(find(ismember(ABDvols(:,1),ABDIr_CellIDs)),:) = ABD.ABDIrContraGradient;
ABD.InterContraGradient(find(ismember(ABDvols(:,1),ABDIc_CellIDs)),:) = ABD.ABDIcContraGradient;
ABD.InterContraGradient(setdiff(1:42,[find(ismember(ABDvols(:,1),ABDIr_CellIDs));find(ismember(ABDvols(:,1),ABDIc_CellIDs))])',:) = [];

ABD.InterRestGradient(find(ismember(ABDvols(:,1),ABDIr_CellIDs)),:) = ABD.ABDIrRestGradient;
ABD.InterRestGradient(find(ismember(ABDvols(:,1),ABDIc_CellIDs)),:) = ABD.ABDIcRestGradient;
ABD.InterRestGradient(setdiff(1:42,[find(ismember(ABDvols(:,1),ABDIr_CellIDs));find(ismember(ABDvols(:,1),ABDIc_CellIDs))])',:) = [];


%%
% plots

figure;
subplot(4,4,1)
plot(0.05:0.1:0.95,nanmean(vertcat(ABD.ABDrSaccadicGradient,ABD.ABDcSaccadicGradient)),'color',lightGreen,'LineWidth',2);
hold on
plot(0.05:0.1:0.95,nanmean(vertcat(ABD.ABDIrSaccadicGradient,ABD.ABDIcSaccadicGradient)),'color',lightMagenta,'LineWidth',2);
box off,
set(gca,'XLim',[0,1],'YLim',[0,0.2]);
axis square;
title('Saccadic');
legend({'ABD','ABDi'},'Location','bestoutside');

subplot(4,4,3)
plot(0.05:0.1:0.95,nanmean(vertcat(ABD.ABDrVestibularGradient,ABD.ABDcVestibularGradient)),'color',lightGreen,'LineWidth',2);
hold on
plot(0.05:0.1:0.95,nanmean(vertcat(ABD.ABDIrVestibularGradient,ABD.ABDIcVestibularGradient)),'color',lightMagenta,'LineWidth',2);
box off,
set(gca,'XLim',[0,1],'YLim',[0,0.2]);
axis square;
title('Vestibular');

subplot(4,4,5)
plot(0.05:0.1:0.95,nanmean(vertcat(ABD.ABDrIntegratorGradient,ABD.ABDcIntegratorGradient)),'color',lightGreen,'LineWidth',2);
hold on
plot(0.05:0.1:0.95,nanmean(vertcat(ABD.ABDIrIntegratorGradient,ABD.ABDIcIntegratorGradient)),'color',lightMagenta,'LineWidth',2);
box off,
set(gca,'XLim',[0,1],'YLim',[0,0.2]);
axis square;
title('Integrator');


subplot(4,4,7)
plot(0.05:0.1:0.95,nanmean(vertcat(ABD.ABDrContraGradient,ABD.ABDcContraGradient)),'color',lightGreen,'LineWidth',2);
hold on
plot(0.05:0.1:0.95,nanmean(vertcat(ABD.ABDIrContraGradient,ABD.ABDIcContraGradient)),'color',lightMagenta,'LineWidth',2);
box off,
set(gca,'XLim',[0,1],'YLim',[0,0.2]);
axis square;
title('Contra');


subplot(4,4,9)
plot(0.05:0.1:0.95,nanmean(vertcat(ABD.ABDrRestGradient,ABD.ABDcRestGradient)),'color',lightGreen,'LineWidth',2);
hold on
plot(0.05:0.1:0.95,nanmean(vertcat(ABD.ABDIrRestGradient,ABD.ABDIcRestGradient)),'color',lightMagenta,'LineWidth',2);
box off,
set(gca,'XLim',[0,1],'YLim',[0,0.2]);
axis square;
title('Rest');


figure
cmapABD = colorcet('D2','N',21);
cmapmotor = flip(cmapABD(1:10,:));
cmapInter = cmapABD(11:21,:);


subplot(4,5,1)
imagesc(ABD.ABDrSaccadicGradient)
set(gca,'XTickLabel',[0.2,0.4,0.6,0.8,1]);
colormap(gca,cmapmotor);


subplot(4,5,2)
imagesc(ABD.ABDrVestibularGradient)
set(gca,'XTickLabel',[0.2,0.4,0.6,0.8,1]);
colormap(gca,cmapmotor);

subplot(4,5,3)
imagesc(ABD.ABDrIntegratorGradient)
set(gca,'XTickLabel',[0.2,0.4,0.6,0.8,1]);
colormap(gca,cmapmotor);

subplot(4,5,4)
imagesc(ABD.ABDrContraGradient)
set(gca,'XTickLabel',[0.2,0.4,0.6,0.8,1]);
colormap(gca,cmapmotor);

subplot(4,5,5)
imagesc(ABD.ABDrRestGradient)
set(gca,'XTickLabel',[0.2,0.4,0.6,0.8,1]);
colormap(gca,cmapmotor);




subplot(4,5,6)
imagesc(ABD.ABDcSaccadicGradient)
set(gca,'XTickLabel',[0.2,0.4,0.6,0.8,1]);
colormap(gca,cmapmotor);

subplot(4,5,7)
imagesc(ABD.ABDcVestibularGradient)
set(gca,'XTickLabel',[0.2,0.4,0.6,0.8,1]);
colormap(gca,cmapmotor);

subplot(4,5,8)
imagesc(ABD.ABDcIntegratorGradient)
set(gca,'XTickLabel',[0.2,0.4,0.6,0.8,1]);
colormap(gca,cmapmotor);

subplot(4,5,9)
imagesc(ABD.ABDcContraGradient)
set(gca,'XTickLabel',[0.2,0.4,0.6,0.8,1]);
colormap(gca,cmapmotor);

subplot(4,5,10)
imagesc(ABD.ABDcRestGradient)
set(gca,'XTickLabel',[0.2,0.4,0.6,0.8,1]);
colormap(gca,cmapmotor);



subplot(4,5,11)
imagesc(ABD.ABDIrSaccadicGradient)
set(gca,'XTickLabel',[0.2,0.4,0.6,0.8,1]);
colormap(gca,cmapInter);

subplot(4,5,12)
imagesc(ABD.ABDIrVestibularGradient)
set(gca,'XTickLabel',[0.2,0.4,0.6,0.8,1]);
colormap(gca,cmapInter);

subplot(4,5,13)
imagesc(ABD.ABDIrIntegratorGradient)
set(gca,'XTickLabel',[0.2,0.4,0.6,0.8,1]);
colormap(gca,cmapInter);

subplot(4,5,14)
imagesc(ABD.ABDIrContraGradient)
set(gca,'XTickLabel',[0.2,0.4,0.6,0.8,1]);
colormap(gca,cmapInter);

subplot(4,5,15)
imagesc(ABD.ABDIrRestGradient)
set(gca,'XTickLabel',[0.2,0.4,0.6,0.8,1]);
colormap(gca,cmapInter);



subplot(4,5,16)
imagesc(ABD.ABDIcSaccadicGradient)
set(gca,'XTickLabel',[0.2,0.4,0.6,0.8,1]);
colormap(gca,cmapInter);

subplot(4,5,17)
imagesc(ABD.ABDIcVestibularGradient)
set(gca,'XTickLabel',[0.2,0.4,0.6,0.8,1]);
colormap(gca,cmapInter);

subplot(4,5,18)
imagesc(ABD.ABDIcIntegratorGradient)
set(gca,'XTickLabel',[0.2,0.4,0.6,0.8,1]);
colormap(gca,cmapInter);

subplot(4,5,19)
imagesc(ABD.ABDIcContraGradient)
set(gca,'XTickLabel',[0.2,0.4,0.6,0.8,1]);
colormap(gca,cmapInter);

subplot(4,5,20)
imagesc(ABD.ABDIcRestGradient)
set(gca,'XTickLabel',[0.2,0.4,0.6,0.8,1]);
colormap(gca,cmapInter);

% plot ABD and ABDi lumped, but sorted by volumes

figure 
subplot(4,4,[1,5])
imagesc(ABD.motorSaccadicGradient);
set(gca,'XTickLabel',[0.2,0.4,0.6,0.8,1],'YTickLabel',[102,115,125,130,140,150]);
title('ADB Saccadic');
colormap(gca,cmapmotor);

subplot(4,4,[2,6])
imagesc(ABD.InterSaccadicGradient);
set(gca,'XTickLabel',[0.2,0.4,0.6,0.8,1],'YTick',[2:4:20],'YTickLabel',[115,130,147,149,169]);
title('ADBi Saccadic');
colormap(gca,cmapInter);

subplot(4,4,[3,7])
imagesc(ABD.motorVestibularGradient);
set(gca,'XTickLabel',[0.2,0.4,0.6,0.8,1],'YTickLabel',[102,115,125,130,140,150]);
title('ADB Vestibular');
colormap(gca,cmapmotor);

subplot(4,4,[4,8])
imagesc(ABD.InterVestibularGradient);
set(gca,'XTickLabel',[0.2,0.4,0.6,0.8,1],'YTick',[2:4:20],'YTickLabel',[115,130,147,149,169]);
title('ADBi Vestibular');
colormap(gca,cmapInter);

subplot(4,4,[9,13])
imagesc(ABD.motorIntegratorGradient);
set(gca,'XTickLabel',[0.2,0.4,0.6,0.8,1],'YTickLabel',[102,115,125,130,140,150]);
title('ADB Integrator');
colormap(gca,cmapmotor);

subplot(4,4,[10,14])
imagesc(ABD.InterIntegratorGradient);
set(gca,'XTickLabel',[0.2,0.4,0.6,0.8,1],'YTick',[2:4:20],'YTickLabel',[115,130,147,149,169]);
title('ADB Integrator');
colormap(gca,cmapInter);

subplot(4,4,[11,15])
imagesc(ABD.motorContraGradient);
set(gca,'XTickLabel',[0.2,0.4,0.6,0.8,1],'YTickLabel',[102,115,125,130,140,150]);
title('ADB Contra');
colormap(gca,cmapmotor);

subplot(4,4,[12,16])
imagesc(ABD.InterContraGradient);
set(gca,'XTickLabel',[0.2,0.4,0.6,0.8,1],'YTick',[2:4:20],'YTickLabel',[115,130,147,149,169]);
title('ADB Contra');
colormap(gca,cmapInter);

%% locate axons in the proxiomal region

ABD.ABDrProximalContraAxons = [];

for i = 1:numel(ABDr_CellIDs)
    temp = find(ABDr(i).PathLength(ABDr(i).isContra)/max(Pvec_tree(ABDr(i).Tree{1}))<0.3);
    ABD.ABDrProximalContraAxons = [ABD.ABDrProximalContraAxons;ABDr(i).Contra(temp)];
    clear temp;
end

ABD.ABDrProximalContraAxons = unique(ABD.ABDrProximalContraAxons);
ABD.ABDrProximalContraAxonsOrigin = getOrigin(ABD.ABDrProximalContraAxons);

ABD.ABDcProximalContraAxons = [];
ABD.ABDcDistalContraAxons = [];
for i = 1:numel(ABDc_CellIDs)
    if ~isempty(ABDc(i).Tree)
        temp = find(ABDc(i).PathLength(ABDc(i).isContra)/max(Pvec_tree(ABDc(i).Tree{1}))<0.3);
        ABD.ABDcProximalContraAxons = [ABD.ABDcProximalContraAxons;ABDc(i).Contra(temp)];
        clear temp;
        temp = find(ABDc(i).PathLength(ABDc(i).isContra)/max(Pvec_tree(ABDc(i).Tree{1}))>0.5);
        ABD.ABDcDistalContraAxons = [ABD.ABDcDistalContraAxons;ABDc(i).Contra(temp)];
        clear temp;
    end
end

ABD.ABDcProximalContraAxons = unique(ABD.ABDcProximalContraAxons);
ABD.ABDcDistalContraAxons = unique(ABD.ABDcDistalContraAxons);
ABD.ABDcProximalContraAxonsOrigin = getOrigin(ABD.ABDcProximalContraAxons);

ABD.ABDrOnlyProximal = setdiff(ABD.ABDrProximalContraAxons,ABD.ABDcProximalContraAxons);
ABD.ABDcOnlyProximal = setdiff(ABD.ABDcProximalContraAxons,ABD.ABDrProximalContraAxons);
ABD.ABDrABDcProximal = intersect(ABD.ABDrProximalContraAxons,ABD.ABDcProximalContraAxons);

%transform_swc_AV(ABD.ABDrOnlyProximal,'k',[],true,false);
%transform_swc_AV(ABD.ABDcOnlyProximal,'r',[],false,false);
%transform_swc_AV(ABD.ABDrABDcProximal,'b',[],false,false);

% Remaining axons

ABD.ABDrProximalRestAxons = [];

for i = 1:numel(ABDr_CellIDs)
    temp = find(ABDr(i).PathLength(ABDr(i).isEverythingElse)/max(Pvec_tree(ABDr(i).Tree{1}))<0.3);
    ABD.ABDrProximalRestAxons = [ABD.ABDrProximalRestAxons;ABDr(i).EverythingElse(temp)];
    clear temp;
end
%ABD.ABDrProximalRestAxons(~ismember(ABD.ABDrProximalRestAxons, ABD.cellID))

ABD.ABDcProximalRestAxons = [];
for i = 1:numel(ABDc_CellIDs)
    if ~isempty(ABDc(i).Tree)
        temp = find(ABDc(i).PathLength(ABDc(i).isEverythingElse)/max(Pvec_tree(ABDc(i).Tree{1}))<0.3);
        ABD.ABDcProximalRestAxons = [ABD.ABDcProximalRestAxons;ABDc(i).EverythingElse(temp)];
        clear temp;
    end
end
