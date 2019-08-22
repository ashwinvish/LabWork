%clear;
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
    df = readtable('/Users/ashwin/Google Drive/Zfish/SynapseDetector/04152019.csv');
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

[~,ABD.sortOrder] = sort(ABD.vol);
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
scatter(ABD.Origin(ABD.ABDcLocs), ABD.ABDcvols,25,lightGreen,'filled','MarkerEdgeColor','k')
hold on;
scatter(ABD.Origin(ABD.ABDrLocs), ABD.ABDrvols,25,lightGreen,'filled','MarkerEdgeColor','k')
scatter(ABD.Origin(ABD.ABDIrLocs), ABD.ABDIrvols,25,lightMagenta,'filled','MarkerEdgeColor','k')
scatter(ABD.Origin(ABD.ABDIcLocs), ABD.ABDIcvols,25,lightMagenta,'filled','MarkerEdgeColor','k')
axis square
xlabel('Medial <---> Lateral');
ylabel('Soma volume \mu^3');

subplot(4,4,2);
scatter(ABD.ABDrvols, ABD.ABDrTotalSynapses,25,lightGreen,'filled','MarkerEdgeColor','k')
hold on;
scatter(ABD.ABDcvols, ABD.ABDcTotalSynapses,25,lightGreen,'filled','MarkerEdgeColor','k')
scatter(ABD.ABDIrvols, ABD.ABDIrTotalSynapses,25,lightMagenta,'filled','MarkerEdgeColor','k')
scatter(ABD.ABDIcvols,ABD.ABDIcTotalSynapses,25,lightMagenta,'filled','MarkerEdgeColor','k')
axis square;
ylabel('Synapses');
xlabel('Soma volume \mu^3');

subplot(4,4,4)
scatter3([ABD.Origin(ABD.ABDrLocs);ABD.Origin(ABD.ABDcLocs)], [ABD.Origin(ABD.ABDrLocs,3);ABD.Origin(ABD.ABDcLocs,3)], [ABD.ABDrvols;ABD.ABDcvols],25,lightGreen,'filled','MarkerEdgeColor','k')
hold on
scatter3([ABD.Origin(ABD.ABDIrLocs);ABD.Origin(ABD.ABDIcLocs)], [ABD.Origin(ABD.ABDIrLocs,3);ABD.Origin(ABD.ABDIcLocs,3)], [ABD.ABDIrvols;ABD.ABDIcvols],25,lightMagenta,'filled','MarkerEdgeColor','k')
axis square;
ylabel('Dorsal <---> Ventral');
xlabel('Medial <---> Lateral');
zlabel('Soma volume \mu^3');
%set(gca, 'YDir','reverse');


subplot(4,4,5)
cmap = colorcet('R3','N',length([ABD.ABDrvols;ABD.ABDcvols]));
scatter([ABD.Origin(ABD.ABDrLocs,1);ABD.Origin(ABD.ABDcLocs,1)]+ [ABD.Origin(ABD.ABDrLocs,3);ABD.Origin(ABD.ABDcLocs,3)], [ABD.ABDrvols;ABD.ABDcvols],...
    25, ABDcolor,'filled','MarkerEdgeColor','k')
hold on;
a = [ABD.Origin(ABD.ABDrLocs);ABD.Origin(ABD.ABDcLocs)]+ [ABD.Origin(ABD.ABDrLocs,3);ABD.Origin(ABD.ABDcLocs,3)];
b = [ABD.ABDrvols;ABD.ABDcvols];
b = b(~isnan(a));
a = a(~isnan(a));
f = showfit(ezfit(a,b,'a*x+b'),'fitcolor',ABDcolor,'dispeqboxmode','off','dispfitlegend','on','corrcoefmode','r2');
legend off;
text(max(a),max(b),sprintf('R^2=%0.2f',f.r^2));
axis square;
xlabel('ML+DV');
ylabel('Soma volume \mu^3');
clear a;
clear b;
offsetAxes(gca);
%set(gca, 'YDir','reverse');


subplot(4,4,6);
scatter(ABD.Origin(ABD.ABDrLocs), ABD.ABDcvols,25,ABDcolor,'filled','MarkerEdgeColor','k');
hold on;
scatter(ABD.Origin(ABD.ABDrLocs), ABD.ABDrvols,25,ABDcolor,'filled','MarkerEdgeColor','k');
a = [ABD.Origin(ABD.ABDrLocs);ABD.Origin(ABD.ABDcLocs)];
b = [ABD.ABDrvols;ABD.ABDcvols];
b = b(~isnan(a));
a = a(~isnan(a));
f = showfit(ezfit(a,b,'a*x+b'),'fitcolor',ABDcolor,'dispeqboxmode','off','dispfitlegend','on','corrcoefmode','r2');
legend off;
text(max(a),max(b),sprintf('R^2=%0.2f',f.r^2));
axis square
xlabel('Medial <---> Lateral');
ylabel('Soma volume \mu^3');


subplot(4,4,7)

temp1 = [ABD.Origin(ABD.ABDrLocs,:);ABD.Origin(ABD.ABDcLocs,:)];
temp2 = [ABD.ABDrvols;ABD.ABDcvols];
temp2 = temp2(~isnan(temp1(:,1)));
temp1 = [temp1(~isnan(temp1(:,1))),temp1(~isnan(temp1(:,2)),2),temp1(~isnan(temp1(:,3)),3)];

[~,ABDmotorSorted] = sort(temp2);

scatter(temp2,temp1(:,1));



f = showfit(ezfit(a,b,'a*x+b'),'fitcolor',ABDcolor,'dispeqboxmode','off','dispfitlegend','on','corrcoefmode','r2');


%%
load ABDr.mat
load ABDc.mat
load ABDIr.mat
load ABDIc.mat
%%

for i = 1:numel(ABDr)
    ABD.ABDrSaccadicSynapses(i) = length(ABDr(i).Saccadic);
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
    ABD.SaccPathLengthABDr{i} = ABDr(i).PathLength(ABDr(i).isSaccadic)./max(Pvec_tree(ABDr(i).Tree{1}));
    ABD.VestPathLengthABDr{i} = ABDr(i).PathLength(ABDr(i).isVestibular)./max(Pvec_tree(ABDr(i).Tree{1}));
    ABD.IntPathLengthABDr{i} = ABDr(i).PathLength(ABDr(i).isIntegrator)./max(Pvec_tree(ABDr(i).Tree{1}));
    ABD.ContraPathLengthABDr{i} = ABDr(i).PathLength(ABDr(i).isContra)./max(Pvec_tree(ABDr(i).Tree{1}));
    ABD.RestPathLengthABDr{i} = ABDr(i).PathLength(ABDr(i).isEverythingElse)./max(Pvec_tree(ABDr(i).Tree{1}));
    
    ABD.ABDrSaccadicGradient(i,:) =  histcounts(ABD.SaccPathLengthABDr{i},0:0.1:1);
    ABD.ABDrVestibularGradient(i,:) =  histcounts(ABD.VestPathLengthABDr{i},0:0.1:1);
    ABD.ABDrIntegratorGradient(i,:) =  histcounts(ABD.IntPathLengthABDr{i},0:0.1:1);
    ABD.ABDrContraGradient(i,:) =  histcounts(ABD.ContraPathLengthABDr{i},0:0.1:1);
    ABD.ABDrRestGradient(i,:) =  histcounts(ABD.RestPathLengthABDr{i},0:0.1:1);   
end
%cellfun(@(x) histogram(x,0:0.1:1),ABD.SaccPathLengthABDr)

% ABDc

for i = 1:numel(ABDc_CellIDs)
    if ~isempty(ABDc(i).Tree)
        ABD.SaccPathLengthABDc{i} = ABDc(i).PathLength(ABDc(i).isSaccadic)./max(Pvec_tree(ABDc(i).Tree{1}));
        ABD.VestPathLengthABDc{i} = ABDc(i).PathLength(ABDc(i).isVestibular)./max(Pvec_tree(ABDc(i).Tree{1}));
        ABD.IntPathLengthABDc{i} = ABDc(i).PathLength(ABDc(i).isIntegrator)./max(Pvec_tree(ABDc(i).Tree{1}));
        ABD.ContraPathLengthABDc{i} = ABDc(i).PathLength(ABDc(i).isContra)./max(Pvec_tree(ABDc(i).Tree{1}));
        ABD.RestPathLengthABDc{i} = ABDc(i).PathLength(ABDc(i).isEverythingElse)./max(Pvec_tree(ABDc(i).Tree{1}));
        
        ABD.ABDcSaccadicGradient(i,:) =  histcounts(ABD.SaccPathLengthABDc{i},0:0.1:1);
        ABD.ABDcVestibularGradient(i,:) =  histcounts(ABD.VestPathLengthABDc{i},0:0.1:1);
        ABD.ABDcIntegratorGradient(i,:) =  histcounts(ABD.IntPathLengthABDc{i},0:0.1:1);
        ABD.ABDcContraGradient(i,:) =  histcounts(ABD.ContraPathLengthABDc{i},0:0.1:1);
        ABD.ABDcRestGradient(i,:) =  histcounts(ABD.RestPathLengthABDc{i},0:0.1:1);
    else
        ABD.ABDcSaccadicGradient(i,:) = repmat(NaN,1,10);
        ABD.ABDcVestibularGradient(i,:) = repmat(NaN,1,10);
        ABD.ABDcIntegratorGradient(i,:) = repmat(NaN,1,10);
        ABD.ABDcContraGradient(i,:) = repmat(NaN,1,10);
        ABD.ABDcRestGradient(i,:)= repmat(NaN,1,10);
    end
end


% combine ABDr and ABDc
ABD.motorCellIDs  = [ABDr_CellIDs';ABDc_CellIDs'];
%ABD.motorCellIDsSorted = intersect(ABD.sortOrder,find(ismember(ABD.cellID,Allmotor)),'stable');

ABD.motorSaccadicGradient = [horzcat(ABD.ABDrvols,ABD.ABDrSaccadicGradient);horzcat(ABD.ABDcvols,ABD.ABDcSaccadicGradient)];
%ABD.motorSaccadicGradient = sortrows(ABD.motorSaccadicGradient,1);
ABD.motorVols = [ABD.vol(ABD.ABDrLocs);ABD.vol(ABD.ABDcLocs)];


ABD.motorVestibularGradient = [horzcat(ABD.ABDrvols,ABD.ABDrVestibularGradient);horzcat(ABD.ABDcvols,ABD.ABDcVestibularGradient)];
%ABD.motorVestibularGradient = sortrows(ABD.motorVestibularGradient,1);

ABD.motorIntegratorGradient = [horzcat(ABD.ABDrvols,ABD.ABDrIntegratorGradient);horzcat(ABD.ABDcvols,ABD.ABDcIntegratorGradient)];
%ABD.motorIntegratorGradient = sortrows(ABD.motorIntegratorGradient,1);

ABD.motorContraGradient = [horzcat(ABD.ABDrvols,ABD.ABDrContraGradient);horzcat(ABD.ABDcvols,ABD.ABDcContraGradient)];
%ABD.motorContraGradient = sortrows(ABD.motorContraGradient,1);

ABD.motorRestGradient = [horzcat(ABD.ABDrvols,ABD.ABDrRestGradient);horzcat(ABD.ABDcvols,ABD.ABDcRestGradient)];
%ABD.motorRestGradient = sortrows(ABD.motorRestGradient,1);
%%
% ABDIr

for i = 1:numel(ABDIr_CellIDs)
    ABD.SaccPathLengthABDIr{i} = ABDIr(i).PathLength(ABDIr(i).isSaccadic)./max(Pvec_tree(ABDIr(i).Tree{1}));
    ABD.VestPathLengthABDIr{i} = ABDIr(i).PathLength(ABDIr(i).isVestibular)./max(Pvec_tree(ABDIr(i).Tree{1}));
    ABD.IntPathLengthABDIr{i} = ABDIr(i).PathLength(ABDIr(i).isIntegrator)./max(Pvec_tree(ABDIr(i).Tree{1}));
    ABD.ContraPathLengthABDIr{i} = ABDIr(i).PathLength(ABDIr(i).isContra)./max(Pvec_tree(ABDIr(i).Tree{1}));
    ABD.RestPathLengthABDIr{i} = ABDIr(i).PathLength(ABDIr(i).isEverythingElse)./max(Pvec_tree(ABDIr(i).Tree{1}));
    
    ABD.ABDIrSaccadicGradient(i,:) =  histcounts(ABD.SaccPathLengthABDIr{i},0:0.1:1);
    ABD.ABDIrVestibularGradient(i,:) =  histcounts(ABD.VestPathLengthABDIr{i},0:0.1:1);
    ABD.ABDIrIntegratorGradient(i,:) =  histcounts(ABD.IntPathLengthABDIr{i},0:0.1:1);
    ABD.ABDIrContraGradient(i,:) =  histcounts(ABD.ContraPathLengthABDIr{i},0:0.1:1);
    ABD.ABDIrRestGradient(i,:) =  histcounts(ABD.RestPathLengthABDIr{i},0:0.1:1);
end


% ABDIc

for i = 1:numel(ABDIc_CellIDs)
    ABD.SaccPathLengthABDIc{i} = ABDIc(i).PathLength(ABDIc(i).isSaccadic)./max(Pvec_tree(ABDIc(i).Tree{1}));
    ABD.VestPathLengthABDIc{i} = ABDIc(i).PathLength(ABDIc(i).isVestibular)./max(Pvec_tree(ABDIc(i).Tree{1}));
    ABD.IntPathLengthABDIc{i} = ABDIc(i).PathLength(ABDIc(i).isIntegrator)./max(Pvec_tree(ABDIc(i).Tree{1}));
    ABD.ContraPathLengthABDIc{i} = ABDIc(i).PathLength(ABDIc(i).isContra)./max(Pvec_tree(ABDIc(i).Tree{1}));
    ABD.RestPathLengthABDIc{i} = ABDIc(i).PathLength(ABDIc(i).isEverythingElse)./max(Pvec_tree(ABDIc(i).Tree{1}));
    
    ABD.ABDIcSaccadicGradient(i,:) =  histcounts(ABD.SaccPathLengthABDIc{i},0:0.1:1);
    ABD.ABDIcVestibularGradient(i,:) =  histcounts(ABD.VestPathLengthABDIc{i},0:0.1:1);
    ABD.ABDIcIntegratorGradient(i,:) =  histcounts(ABD.IntPathLengthABDIc{i},0:0.1:1);
    ABD.ABDIcContraGradient(i,:) =  histcounts(ABD.ContraPathLengthABDIc{i},0:0.1:1);
    ABD.ABDIcRestGradient(i,:) =  histcounts(ABD.RestPathLengthABDIc{i},0:0.1:1);
end




% combine ABDIr and ABDIc
ABD.InterCellIDs = [ABDIr_CellIDs';ABDIc_CellIDs'];
%ABD.InterCellIDsSorted = intersect(ABD.sortOrder,find(setdiff(ABD.cellID,Allmotor)),'stable');

ABD.InterSaccadicGradient = [horzcat(ABD.ABDIrvols,ABD.ABDIrSaccadicGradient);horzcat(ABD.ABDIcvols,ABD.ABDIcSaccadicGradient)];
ABD.InterSaccadicGradient = sortrows(ABD.InterSaccadicGradient,1);

ABD.InterVestibularGradient = [horzcat(ABD.ABDIrvols,ABD.ABDIrVestibularGradient);horzcat(ABD.ABDIcvols,ABD.ABDIcVestibularGradient)];
ABD.InterVestibularGradient = sortrows(ABD.InterVestibularGradient,1);

ABD.InterIntegratorGradient = [horzcat(ABD.ABDIrvols,ABD.ABDIrIntegratorGradient);horzcat(ABD.ABDIcvols,ABD.ABDIcIntegratorGradient)];
ABD.InterIntegratorGradient = sortrows(ABD.InterIntegratorGradient,1);

ABD.InterContraGradient = [horzcat(ABD.ABDIrvols,ABD.ABDIrContraGradient);horzcat(ABD.ABDIcvols,ABD.ABDIcContraGradient)];
ABD.InterContraGradient = sortrows(ABD.InterContraGradient,1);

ABD.InterRestGradient = [horzcat(ABD.ABDIrvols,ABD.ABDIrRestGradient);horzcat(ABD.ABDIcvols,ABD.ABDIcRestGradient)];
ABD.InterRestGradient = sortrows(ABD.InterRestGradient ,1);

%%
figure;


subplot(4,4,1)
shadedErrorBar(0.1:0.1:1,nanmean(ABD.motorSaccadicGradient(:,2:end)),nanstd(ABD.motorSaccadicGradient(:,2:end))./sqrt(29),'lineprops',{'Color',ABDcolor,'Linewidth',2});
hold on;
shadedErrorBar(0.1:0.1:1,nanmean(ABD.InterSaccadicGradient(:,2:end)),nanstd(ABD.InterSaccadicGradient(:,2:end))./sqrt(21),'lineprops',{'Color',ABDicolor,'Linewidth',2});
axis square;
title('Saccadic axons');
legend({'M','I'},'Location','bestoutside');
xlabel('Norm. pathlength');
ylabel('Number of synapses');
offsetAxes(gca);

subplot(4,4,2)
shadedErrorBar(0.1:0.1:1,nanmean(ABD.motorVestibularGradient(:,2:end)),nanstd(ABD.motorVestibularGradient(:,2:end))./sqrt(29),'lineprops',{'Color',ABDcolor,'Linewidth',2});
hold on;
shadedErrorBar(0.1:0.1:1,nanmean(ABD.InterVestibularGradient(:,2:end)),nanstd(ABD.InterVestibularGradient(:,2:end))./sqrt(21),'lineprops',{'Color',ABDicolor,'Linewidth',2});
axis square;
title('Vestibular axons');
xlabel('Norm. pathlength');
ylabel('Number of synapses');
legend({'M','I'},'Location','bestoutside');
offsetAxes(gca);

figure;
for i = 1:numel(ABD.motorCellIDs)
subplot(5,6,i)
plot(0.1:0.1:1,ABD.motorSaccadicGradient(i,2:end)./max(ABD.motorSaccadicGradient(i,2:end)));
end


%%
% plots

figure;
subplot(4,4,1)
plot(0.05:0.1:0.95,nanmean(vertcat(ABD.ABDrSaccadicGradient,ABD.ABDcSaccadicGradient)),'color',lightGreen,'LineWidth',2);
hold on
plot(0.05:0.1:0.95,nanmean(vertcat(ABD.ABDIrSaccadicGradient,ABD.ABDIcSaccadicGradient)),'color',lightMagenta,'LineWidth',2);
box off,
%set(gca,'XLim',[0,1],'YLim',[0,0.2]);
axis square;
title('Saccadic');
legend({'ABD','ABDi'},'Location','bestoutside');

subplot(4,4,3)
plot(0.05:0.1:0.95,nanmean(vertcat(ABD.ABDrVestibularGradient,ABD.ABDcVestibularGradient)),'color',lightGreen,'LineWidth',2);
hold on
plot(0.05:0.1:0.95,nanmean(vertcat(ABD.ABDIrVestibularGradient,ABD.ABDIcVestibularGradient)),'color',lightMagenta,'LineWidth',2);
box off,
%set(gca,'XLim',[0,1],'YLim',[0,0.2]);
axis square;
title('Vestibular');

subplot(4,4,5)
plot(0.05:0.1:0.95,nanmean(vertcat(ABD.ABDrIntegratorGradient,ABD.ABDcIntegratorGradient)),'color',lightGreen,'LineWidth',2);
hold on
plot(0.05:0.1:0.95,nanmean(vertcat(ABD.ABDIrIntegratorGradient,ABD.ABDIcIntegratorGradient)),'color',lightMagenta,'LineWidth',2);
box off,
%set(gca,'XLim',[0,1],'YLim',[0,0.2]);
axis square;
title('Integrator');


subplot(4,4,7)
plot(0.05:0.1:0.95,nanmean(vertcat(ABD.ABDrContraGradient,ABD.ABDcContraGradient)),'color',lightGreen,'LineWidth',2);
hold on
plot(0.05:0.1:0.95,nanmean(vertcat(ABD.ABDIrContraGradient,ABD.ABDIcContraGradient)),'color',lightMagenta,'LineWidth',2);
box off,
%set(gca,'XLim',[0,1],'YLim',[0,0.2]);
axis square;
title('Contra');


subplot(4,4,9)
plot(0.05:0.1:0.95,nanmean(vertcat(ABD.ABDrRestGradient,ABD.ABDcRestGradient)),'color',lightGreen,'LineWidth',2);
hold on
plot(0.05:0.1:0.95,nanmean(vertcat(ABD.ABDIrRestGradient,ABD.ABDIcRestGradient)),'color',lightMagenta,'LineWidth',2);
box off,
%set(gca,'XLim',[0,1],'YLim',[0,0.2]);
axis square;
title('Rest');


figure;
cmapABD = colorcet('D2','N',21);
cmapmotor = flip(cmapABD(1:11,:));
cmapInter = cmapABD(11:end,:);


subplot(4,5,1)
heatmap(ABD.ABDrSaccadicGradient,'ColorScaling','scaledrows','CellLabelColor','none',...
   'Colormap',cmapmotor,'YDisplayLabels',round(ABD.ABDrvols),'XDisplayLabels',0.1:0.1:1);


subplot(4,5,2)
heatmap(ABD.ABDrVestibularGradient,'ColorScaling','scaledrows','CellLabelColor',...
    'none','Colormap',cmapmotor,'YDisplayLabels',round(ABD.ABDrvols),'XDisplayLabels',0.1:0.1:1);


subplot(4,5,3)
heatmap(ABD.ABDrIntegratorGradient,'ColorScaling','scaledrows','CellLabelColor',...
    'none','Colormap',cmapmotor,'YDisplayLabels',round(ABD.ABDrvols),'XDisplayLabels',0.1:0.1:1);

subplot(4,5,4)
heatmap(ABD.ABDrContraGradient,'ColorScaling','scaledrows','CellLabelColor',...
    'none','Colormap',cmapmotor,'YDisplayLabels',round(ABD.ABDrvols),'XDisplayLabels',0.1:0.1:1);


subplot(4,5,5)
heatmap(ABD.ABDrRestGradient,'ColorScaling','scaledrows','CellLabelColor',...
    'none','Colormap',cmapmotor,'YDisplayLabels',round(ABD.ABDrvols),'XDisplayLabels',0.1:0.1:1);


subplot(4,5,6)
heatmap(ABD.ABDcSaccadicGradient,'ColorScaling','scaledrows','CellLabelColor',...
    'none','Colormap',cmapmotor,'YDisplayLabels',round(ABD.ABDcvols),'XDisplayLabels',0.1:0.1:1);


subplot(4,5,7)
heatmap(ABD.ABDcVestibularGradient,'ColorScaling','scaledrows','CellLabelColor',...
    'none','Colormap',cmapmotor,'YDisplayLabels',round(ABD.ABDcvols),'XDisplayLabels',0.1:0.1:1);

subplot(4,5,8)
heatmap(ABD.ABDcIntegratorGradient,'ColorScaling','scaledrows','CellLabelColor',...
    'none','Colormap',cmapmotor,'YDisplayLabels',round(ABD.ABDcvols),'XDisplayLabels',0.1:0.1:1);

subplot(4,5,9)
heatmap(ABD.ABDcContraGradient,'ColorScaling','scaledrows','CellLabelColor',...
    'none','Colormap',cmapmotor,'YDisplayLabels',round(ABD.ABDcvols),'XDisplayLabels',0.1:0.1:1);

subplot(4,5,10)
heatmap(ABD.ABDcRestGradient,'ColorScaling','scaledrows','CellLabelColor',...
    'none','Colormap',cmapmotor,'YDisplayLabels',round(ABD.ABDcvols),'XDisplayLabels',0.1:0.1:1);



subplot(4,5,11)
heatmap(ABD.ABDIrSaccadicGradient,'ColorScaling','scaledrows','CellLabelColor',...
    'none','Colormap',cmapInter,'YDisplayLabels',round(ABD.ABDIrvols),'XDisplayLabels',0.1:0.1:1);


subplot(4,5,12)
heatmap(ABD.ABDIrVestibularGradient,'ColorScaling','scaledrows','CellLabelColor',...
    'none','Colormap',cmapInter,'YDisplayLabels',round(ABD.ABDIrvols),'XDisplayLabels',0.1:0.1:1);


subplot(4,5,13)
heatmap(ABD.ABDIrIntegratorGradient,'ColorScaling','scaledrows','CellLabelColor',...
    'none','Colormap',cmapInter,'YDisplayLabels',round(ABD.ABDIrvols),'XDisplayLabels',0.1:0.1:1);


subplot(4,5,14)
heatmap(ABD.ABDIrContraGradient,'ColorScaling','scaledrows','CellLabelColor',...
    'none','Colormap',cmapInter,'YDisplayLabels',round(ABD.ABDIrvols),'XDisplayLabels',0.1:0.1:1);



subplot(4,5,15)
heatmap(ABD.ABDIrRestGradient,'ColorScaling','scaledrows','CellLabelColor',...
    'none','Colormap',cmapInter,'YDisplayLabels',round(ABD.ABDIrvols),'XDisplayLabels',0.1:0.1:1);



subplot(4,5,16)
heatmap(ABD.ABDIcSaccadicGradient,'ColorScaling','scaledrows','CellLabelColor',...
    'none','Colormap',cmapInter,'YDisplayLabels',round(ABD.ABDIcvols),'XDisplayLabels',0.1:0.1:1);


subplot(4,5,17)
heatmap(ABD.ABDIcVestibularGradient,'ColorScaling','scaledrows','CellLabelColor',...
    'none','Colormap',cmapInter,'YDisplayLabels',round(ABD.ABDIcvols),'XDisplayLabels',0.1:0.1:1);


subplot(4,5,18)
heatmap(ABD.ABDIcIntegratorGradient,'ColorScaling','scaledrows','CellLabelColor',...
    'none','Colormap',cmapInter,'YDisplayLabels',round(ABD.ABDIcvols),'YDisplayLabels',0.1:0.1:1);


subplot(4,5,19)
heatmap(ABD.ABDIcContraGradient,'ColorScaling','scaledrows','CellLabelColor',...
    'none','Colormap',cmapInter,'YDisplayLabels',round(ABD.ABDIcvols),'XDisplayLabels',0.1:0.1:1);



subplot(4,5,20)
heatmap(ABD.ABDIcRestGradient,'ColorScaling','scaledrows','CellLabelColor',...
    'none','Colormap',cmapInter,'YDisplayLabels',round(ABD.ABDIcvols),'XDisplayLabels',0.1:0.1:1);


% plot ABD and ABDi lumped, but sorted by volumes

figure;
ABD.motorVols = ABD.motorSaccadicGradient(:,1);

subplot(4,4,1)
heatmap(ABD.motorSaccadicGradient(:,2:end),'ColorScaling','scaledrows','CellLabelColor',...
    'none','Colormap',cmapmotor,'YDisplayLabels',round(ABD.motorVols),'XDisplayLabels',0.1:0.1:1);
title('ADB Saccadic');


% hold on;
% [~,b] = max(ABD.motorSaccadicGradient(:,2:end),[],2);
% plot(b,'ko');
% showfit(ezfit(1:size(b,1),b,'a*x+b'),'fitcolor','k','dispeqboxmode','off','dispfitlegend','on','corrcoefmode','r2');
% clear b;

ABD.interVols = ABD.InterSaccadicGradient(:,1);
subplot(4,4,2)
heatmap(ABD.InterSaccadicGradient(:,2:end),'ColorScaling','scaledrows','CellLabelColor',...
    'none','Colormap',cmapInter,'YDisplayLabels',round(ABD.interVols),'XDisplayLabels',0.1:0.1:1);
title('ADBi Saccadic');
%daspect([2,1,1]);
colormap(gca,cmapInter);
% hold on;
% [~,b] = max(ABD.InterSaccadicGradient(:,2:end),[],2);
% plot(11-b,'ko');
% showfit(ezfit(1:size(b,1),11-b,'a*x'),'fitcolor','k','dispeqboxmode','off','corrcoefmode','r2');
% clear b;

subplot(4,4,3)
heatmap(ABD.motorVestibularGradient(:,2:end),'ColorScaling','scaledrows','CellLabelColor',...
    'none','Colormap',cmapmotor,'YDisplayLabels',round(ABD.motorVols),'XDisplayLabels',0.1:0.1:1);

title('ADB Vestibular');
colormap(gca,cmapmotor);

subplot(4,4,4)
heatmap(ABD.InterVestibularGradient(:,2:end),'ColorScaling','scaledrows','CellLabelColor',...
    'none','Colormap',cmapInter,'YDisplayLabels',round(ABD.interVols),'XDisplayLabels',0.1:0.1:1);
title('ADBi Vestibular');
colormap(gca,cmapInter);

subplot(4,4,5)
heatmap(ABD.motorIntegratorGradient(:,2:end),'ColorScaling','scaledrows','CellLabelColor',...
    'none','Colormap',cmapmotor,'YDisplayLabels',round(ABD.motorVols),'XDisplayLabels',0.1:0.1:1);
title('ADB Integrator');
colormap(gca,cmapmotor);

subplot(4,4,6)
heatmap(ABD.InterIntegratorGradient(:,2:end),'ColorScaling','scaledrows','CellLabelColor',...
    'none','Colormap',cmapInter,'YDisplayLabels',round(ABD.interVols),'XDisplayLabels',0.1:0.1:1);
title('ADB Integrator');
colormap(gca,cmapInter);

subplot(4,4,7)
heatmap(ABD.motorContraGradient(:,2:end),'ColorScaling','scaledrows','CellLabelColor',...
    'none','Colormap',cmapmotor,'YDisplayLabels',round(ABD.motorVols),'XDisplayLabels',0.1:0.1:1);
title('ADB Contra');
colormap(gca,cmapmotor);

subplot(4,4,8)
heatmap(ABD.InterContraGradient(:,2:end),'ColorScaling','scaledcolumns','CellLabelColor',...
    'none','Colormap',cmapInter,'YDisplayLabels',round(ABD.interVols),'XDisplayLabels',0.1:0.1:1);
title('ADBi Contra');
colormap(gca,cmapInter);

%% Arrange Medial to Lateral

ABD.MotorOrigins = [vertcat(ABDr.Origin);vertcat(ABDc.Origin)];
ABD.motorVols(15) = [];
ABD.motorVols(29) = [];
ABD.motorSaccadicGradient(15,:) = [];
ABD.motorSaccadicGradient(29,:) = [];

figure;
% arranged Small to Big
 [~,ABDMotorSmallBig] = sort(ABD.motorVols);

subplot(4,4,[1,5])
heatmap(ABD.motorSaccadicGradient(ABDMotorSmallBig,2:end),'ColorScaling','scaledrows','CellLabelColor',...
    'none','Colormap',cmapmotor,'YDisplayLabels',round(ABD.motorVols(ABDMotorSmallBig)),'XDisplayLabels',0.1:0.1:1);
title('Small to Big');

% arranged ML+DV
[~,ABDMotorMLDVSort] = sort(ABD.MotorOrigins(:,1) + ABD.MotorOrigins(:,3));

subplot(4,4,[2,6])
heatmap(ABD.motorSaccadicGradient(ABDMotorMLDVSort,2:end),'ColorScaling','scaledrows','CellLabelColor',...
    'none','Colormap',cmapmotor,'YDisplayLabels',round(ABD.motorVols(ABDMotorMLDVSort)),'XDisplayLabels',0.1:0.1:1);
title('ML+DV');


subplot(4,4,9)
ToptwothirdMean = nanmean(ABD.motorSaccadicGradient(ABDMotorMLDVSort(1:10),2:end));
ToptwothirdSEM = nanstd(ABD.motorSaccadicGradient(ABDMotorMLDVSort(1:10),2:end))./ sqrt(9);

LastonethirdMean = nanmean(ABD.motorSaccadicGradient(ABDMotorMLDVSort(22:end),2:end));
LastonethirdSEM = nanstd(ABD.motorSaccadicGradient(ABDMotorMLDVSort(22:end),2:end))./sqrt(9);


shadedErrorBar(0.1:0.1:1,ToptwothirdMean,ToptwothirdSEM,'lineprops',{'Color','r'});
hold on;
shadedErrorBar(0.1:0.1:1,LastonethirdMean,LastonethirdSEM,'lineprops',{'Color','k'});
ylabel('average synapses');
legend({'top 1/3','bottom 1/3'},'location','best');


% un normalized
subplot(4,4,13)
ToptwothirdMean = nanmean(ABD.motorSaccadicGradient(ABDMotorMLDVSort(1:21),2:end));
ToptwothirdSEM = nanstd(ABD.motorSaccadicGradient(ABDMotorMLDVSort(1:21),2:end))./ sqrt(20);

LastonethirdMean = nanmean(ABD.motorSaccadicGradient(ABDMotorMLDVSort(22:end),2:end));
LastonethirdSEM = nanstd(ABD.motorSaccadicGradient(ABDMotorMLDVSort(22:end),2:end))./sqrt(9);

shadedErrorBar(0.1:0.1:1,ToptwothirdMean,ToptwothirdSEM,'lineprops',{'Color','r'});
hold on;
shadedErrorBar(0.1:0.1:1,LastonethirdMean,LastonethirdSEM,'lineprops',{'Color','k'});
xlabel('Norm. path length');
legend({'top 2/3','bottom 1/3'},'location','best');


% vestibular nuerons
ABD.motorVestibularGradient(15,:) = [];
ABD.motorVestibularGradient(29,:) = [];

subplot(4,4,[2,6])
heatmap(ABD.motorVestibularGradient(ABDMotorMLDVSort,2:end),'ColorScaling','scaledrows','CellLabelColor',...
    'none','Colormap',cmapmotor,'YDisplayLabels',round(ABD.motorVols(ABDMotorMLDVSort)),'XDisplayLabels',0.1:0.1:1);
title('ML+DV');


subplot(4,4,14)
ToptwothirdMean = nanmean(ABD.motorVestibularGradient(ABDMotorMLDVSort(1:21),2:end));
ToptwothirdSEM = nanstd(ABD.motorVestibularGradient(ABDMotorMLDVSort(1:21),2:end))./ sqrt(20);

LastonethirdMean = nanmean(ABD.motorVestibularGradient(ABDMotorMLDVSort(22:end),2:end));
LastonethirdSEM = nanstd(ABD.motorVestibularGradient(ABDMotorMLDVSort(22:end),2:end))./sqrt(9);

shadedErrorBar(0.1:0.1:1,ToptwothirdMean,ToptwothirdSEM,'lineprops',{'Color','r'});
hold on;
shadedErrorBar(0.1:0.1:1,LastonethirdMean,LastonethirdSEM,'lineprops',{'Color','k'});



% plot for R and C seperately

% [~,ABD.ABDrMLDVorder] = sort(ABD.MotorOrigins(1:14,1)+ABD.MotorOrigins(1:14,3));
% 
% 
% subplot(4,4,[2 6])
% heatmap(ABD.motorSaccadicGradient(ABD.ABDrMLDVorder,2:end),'ColorScaling','scaledrows','CellLabelColor',...
%     'none','Colormap',cmapmotor,'YDisplayLabels',round(ABD.motorVols(ABD.ABDrMLDVorder)),'XDisplayLabels',0.1:0.1:1);
% title('ABDr: ML+DV');
% 
% subplot(4,4,10)
% % plot(sort(ABD.MotorOrigins(:,1) + ABD.MotorOrigins(:,3)),round(ABD.motorVols(ABDMotorMLDVSort)),'o', ...
% %     'MarkerFaceColor',lightGreen,'MarkerEdgeColor','none');
%  axis square;
% % box off;
% % xlabel('ML+DV');
% % ylabel('Soma volume');
% 
% ToptwothirdMean = nanmean(ABD.motorSaccadicGradient(ABD.ABDrMLDVorder(1:9),2:end));
% ToptwothirdSEM = nanstd(ABD.motorSaccadicGradient(ABD.ABDrMLDVorder(1:9),2:end))./ sqrt(8);
% 
% LastonethirdMean = nanmean(ABD.motorSaccadicGradient(ABD.ABDrMLDVorder(10:end),2:end));
% LastonethirdSEM = nanstd(ABD.motorSaccadicGradient(ABD.ABDrMLDVorder(10:end),2:end))./sqrt(5);
% 
% subplot(4,4,14)
% shadedErrorBar(0.1:0.1:1,ToptwothirdMean,ToptwothirdSEM,'lineprops',{'Color','r'});
% hold on;
% shadedErrorBar(0.1:0.1:1,LastonethirdMean,LastonethirdSEM,'lineprops',{'Color','k'});
% 
% 
% 
% [~,ABD.ABDcMLDVorder] = sort(ABD.MotorOrigins(15:end,1)+ABD.MotorOrigins(15:end,3));
% 
% subplot(4,4,[3 7])
% heatmap(ABD.motorSaccadicGradient(14+ABD.ABDcMLDVorder,2:end),'ColorScaling','scaledrows','CellLabelColor',...
%     'none','Colormap',cmapmotor,'YDisplayLabels',round(ABD.motorVols(14+ABD.ABDcMLDVorder)),'XDisplayLabels',0.1:0.1:1);
% title('ABDc: ML+DV');
% 
% subplot(4,4,11)
% % plot(sort(ABD.MotorOrigins(:,1) + ABD.MotorOrigins(:,3)),round(ABD.motorVols(ABDMotorMLDVSort)),'o', ...
% %     'MarkerFaceColor',lightGreen,'MarkerEdgeColor','none');
%  axis square;
% % box off;
% % xlabel('ML+DV');
% % ylabel('Soma volume');
% 
% ToptwothirdMean = nanmean(ABD.motorSaccadicGradient(14+ABD.ABDcMLDVorder(1:9),2:end));
% ToptwothirdSEM = nanstd(ABD.motorSaccadicGradient(14+ABD.ABDcMLDVorder(1:9),2:end))./ sqrt(9);
% 
% LastonethirdMean = nanmean(ABD.motorSaccadicGradient(14+ABD.ABDcMLDVorder(10:end),2:end));
% LastonethirdSEM = nanstd(ABD.motorSaccadicGradient(14+ABD.ABDcMLDVorder(10:end),2:end))./sqrt(7);
% 
% subplot(4,4,15)
% shadedErrorBar(0.1:0.1:1,ToptwothirdMean,ToptwothirdSEM,'lineprops',{'Color','r'});
% hold on;
% shadedErrorBar(0.1:0.1:1,LastonethirdMean,LastonethirdSEM,'lineprops',{'Color','k'});
% 
% 
% 
% 

%% Plot Neuros with ML+DV order

cmap = colorcet('R3','N',numel(ABD.motorCellIDs));
transform_swc_AV(ABD.motorCellIDs(ABDMotorMLDVSort), cmap,[],true,false);
colormap(cmap);
colorbar('Ticks',[0,1],'TickLabels',{'Medio-dorsal','Ventro-lateral'});

%% sort by mono Bi order

ABDc_monostratified = [82212,77654,77295,77646]; % removed 82213
ABDc_bistratified = [77662,77154,77688,77292,77658,77657,77682,77296,82195,77652,82197,81172,77628,82196];

ABDr_monostratified = [82145,77302,82140,82143,82146]; % removed 77299
ABDr_bistratified = [82192,77709,77672,82193,77305,77710,77300,77648,82194,77301,77705];

ABDStratCellIDs = [ABDr_monostratified'; ABDc_monostratified';ABDr_bistratified'; ABDc_bistratified'];
[~,ABDStratOrder] = ismember(ABDStratCellIDs,ABD.motorCellIDs);

ABDStratOrder = ABDStratOrder(ABDStratOrder~=0); 

ABD.motorSaccadicGradientNormalized = ABD.motorSaccadicGradient(:,2:end)./max(ABD.motorSaccadicGradient(:,2:end),[],2);

subplot(4,4,[1,5])
heatmap(ABD.motorSaccadicGradient(:,2:end),'ColorScaling','scaledrows','CellLabelColor',...
    'none','Colormap',cmapmotor,'YDisplayLabels',round(ABD.motorVols),'XDisplayLabels',0.1:0.1:1);
title('size');

subplot(4,4,9)
Quantile25th = find(ABD.motorSaccadicGradient(:,1)<quantile(ABD.motorSaccadicGradient(:,1),0.25));
Quantile25thRest = find(ABD.motorSaccadicGradient(:,1)>quantile(ABD.motorSaccadicGradient(:,1),0.25));

Quantile25thMean = nanmean(ABD.motorSaccadicGradient(Quantile25th,2:end));
Quantile25thSEM = nanstd(ABD.motorSaccadicGradient(Quantile25th,2:end))./sqrt(length(Quantile25th));

Quantile25thRestMean = nanmean(ABD.motorSaccadicGradient(Quantile25thRest,2:end));
Quantile25thRestSEM = nanstd(ABD.motorSaccadicGradient(Quantile25thRest,2:end))./sqrt(length(Quantile25thRest));

shadedErrorBar(0.1:0.1:1,Quantile25thMean,Quantile25thSEM,'lineprops',{'Color','k'});
hold on;
shadedErrorBar(0.1:0.1:1,Quantile25thRestMean,Quantile25thRestSEM,'lineprops',{'Color','r'});
set(gca, 'XLim',[0.1,1]);
title('25th Quantile and rest');

subplot(4,4,13)

Quantile75th = find(ABD.motorSaccadicGradient(:,1)>quantile(ABD.motorSaccadicGradient(:,1),0.75));
Quantile75thRest = find(ABD.motorSaccadicGradient(:,1)<quantile(ABD.motorSaccadicGradient(:,1),0.75));

Quantile75thMean = nanmean(ABD.motorSaccadicGradient(Quantile75th,2:end));
Quantile75thSEM = nanstd(ABD.motorSaccadicGradient(Quantile75th,2:end))./sqrt(length(Quantile75th));

Quantile75thRestMean = nanmean(ABD.motorSaccadicGradient(Quantile75thRest,2:end));
Quantile75thRestSEM = nanstd(ABD.motorSaccadicGradient(Quantile75thRest,2:end)) ./ sqrt(length(Quantile75thRest));

shadedErrorBar(0.1:0.1:1,Quantile75thMean,Quantile75thSEM,'lineprops',{'Color','k'});
hold on;
shadedErrorBar(0.1:0.1:1,Quantile75thRestMean,Quantile75thRestSEM,'lineprops',{'Color','r'});
set(gca, 'XLim',[0.1,1]);
title('75th Quantile and rest');


%offsetAxes(gca);

subplot(4,4,[3,7])
heatmap(ABD.motorSaccadicGradient(ABDStratOrder,2:end),'ColorScaling','scaledrows','CellLabelColor',...
    'none','Colormap',cmapmotor,'YDisplayLabels',round(ABD.motorVols(ABDStratOrder)),'XDisplayLabels', 0.1:0.1:1 );
title('mono_bi');

monoMean = nanmean(ABD.motorSaccadicGradient(ABDStratOrder(2:9),2:end));
monoSEM = nanstd(ABD.motorSaccadicGradient(ABDStratOrder(1:9),2:end))./sqrt(8);

BiMean = nanmean(ABD.motorSaccadicGradient(ABDStratOrder(10:end),2:end));
BiSEM = nanstd(ABD.motorSaccadicGradient(ABDStratOrder(10:end),2:end))./sqrt(21);

subplot(4,4,11)
shadedErrorBar(0.1:0.1:1,monoMean,monoSEM,'lineprops',{'Color','k'});
hold on;
shadedErrorBar(0.1:0.1:1,BiMean,BiSEM,'lineprops',{'Color','r'});
set(gca, 'XLim',[0.1,1]);


ABDcMyelinated = [77292,77296,77154,77295,81172];


%% translate ABD neurons to same position to see effect

% get Motor neuron Origins



for i = 1:numel(ABD.motorCellIDs)
    if ismember(ABD.motorCellIDs(i), ABDr_CellIDs)
        [~,b] = ismember(ABD.motorCellIDs(i), ABDr_CellIDs);
        ABD.MotorSaccadicSites{i} = ABDr(b).PreSynCoordsTransformed(ABDr(b).isSaccadic,:);
    elseif ismember(ABD.motorCellIDs(i), ABDc_CellIDs)
        [~,b] = ismember(ABD.motorCellIDs(i), ABDc_CellIDs);
        ABD.MotorSaccadicSites{i} = ABDc(b).PreSynCoordsTransformed(ABDc(b).isSaccadic,:);
    end
end


newOrigin = ABD.MotorOrigins(ABDMotorMLSort(1),:);


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

%% plot numbers for all neurons





