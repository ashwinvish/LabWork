clear;
addpath(genpath('/Users/ashwin/Documents/'));

load ABDPutativeSaccadic.mat
load ABDiPutativeSaccadic.mat
load ABD_ABDi_PutativeSaccadic.mat

lateralVSaccadic = [80163 80167 80177 80179 76688 80204 80206 80210 76682];



colorPallete = cbrewer('div','BrBG',5);
colorSchemes

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
scatter(ABD.Origin(ABD.ABDrLocs), ABD.ABDrvols,25,ABDcolor,'filled','MarkerEdgeColor','k');
hold on;
scatter(ABD.Origin(ABD.ABDcLocs), ABD.ABDcvols,25,ABDcolor,'filled','MarkerEdgeColor','k');
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
clear a;
clear b;


subplot(4,4,7)
a = [ABD.Origin(ABD.ABDrLocs);ABD.Origin(ABD.ABDcLocs)]+ [ABD.Origin(ABD.ABDrLocs,3);ABD.Origin(ABD.ABDcLocs,3)];
histogram(a,min(a):6:max(a),'FaceColor',ABDcolor,'EdgeColor','none');
hold on;
xlabel({'Dorsomedial <--> Ventrolateral',' ML+DV'});
ylabel('count')
yyaxis('right');
scatter([ABD.Origin(ABD.ABDrLocs,1);ABD.Origin(ABD.ABDcLocs,1)]+ [ABD.Origin(ABD.ABDrLocs,3);ABD.Origin(ABD.ABDcLocs,3)], [ABD.ABDrvols;ABD.ABDcvols],...
    25, ABDcolor,'filled','MarkerEdgeColor','k');
set(gca,'YLim',[75,160],'YColor','k');
% histogram(a,min(a):6:max(a),'Normalization','cumcount','DisplayStyle','stairs','LineWidth',2,'EdgeColor','r');
line([437,437],[0,160],'color','k','LineStyle',':');
 ylabel('Soma volume \mu^3');
box off;
axis square;
clear a;

subplot(4,4,8)
scatter(ABD.ABDrvols, ABD.ABDrTotalSynapses,25,ABDcolor,'filled','MarkerEdgeColor',ABDcolor)
hold on;
scatter(ABD.ABDcvols, ABD.ABDcTotalSynapses,25,lightGreen,'MarkerEdgeColor',lightGreen)
ylabel('Synapses');
xlabel('Soma volume \mu^3');
axis square;
box off;
legend({'r5','r6'},'Location','bestoutside');


% subplot(4,4,8)
% map = cbrewer('seq','Greens',14); % cell are already sorted accordig to size
% scatter([ABD.Origin(ABD.ABDrLocs,1)+ABD.Origin(ABD.ABDrLocs,3)], ABD.ABDrTotalSynapses,50,map(sortr5,:),'filled','MarkerEdgeColor','k')
% ylabel('Synapses');
% xlabel('ML+DV ');
% %axis square;
% daspect([1,5,1])
% colormap(map)
% colorbar;
% caxis([ABD.ABDrvols(1), ABD.ABDrvols(end)]);
% legend('r5','Location','bestoutside');

subplot(4,4,9)
map = cbrewer('seq','Greens',17);
scatter([ABD.Origin(ABD.ABDcLocs,1)+ABD.Origin(ABD.ABDcLocs,3)], ABD.ABDcTotalSynapses,50,map,'filled','MarkerEdgeColor','k')
ylabel('Synapses');
xlabel('ML+DV');
daspect([1,5,1])
box off;
colormap(map)
colorbar;
caxis([ABD.ABDcvols(1), ABD.ABDcvols(end)]);
legend('r6','Location','bestoutside');

figure;
subplot(4,4,1)
map = cbrewer('seq','Greens',31);
[~,sortr5r6] = sort([ABD.ABDrvols;ABD.ABDcvols]) ;
scatter([ABD.Origin(ABD.ABDrLocs,1)+ABD.Origin(ABD.ABDrLocs,3); ABD.Origin(ABD.ABDcLocs,1)+ABD.Origin(ABD.ABDcLocs,3)],...
   [ABD.ABDrTotalSynapses; ABD.ABDcTotalSynapses],50,map(sortr5r6,:),'filled','MarkerEdgeColor','k')
ylabel('Synapses');
xlabel('ML+DV');
axis square;
box off;
colormap(map)
colorbar;
caxis([ABD.ABDcvols(1), ABD.ABDcvols(end)]);

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

ABD.motorSaccadicSynapseNumbers = [ABD.ABDrSaccadicSynapses,ABD.ABDcSaccadicSynapses];
ABD.motorVestibularSynapseNumbers = [ABD.ABDrVestibularSynapses,ABD.ABDcVestibularSynapses];
ABD.motorIntegratorSynapseNumbers = [ABD.ABDrIntegratorSynapses,ABD.ABDcIntegratorSynapses];
ABD.motorContraSynapseNumbers = [ABD.ABDrContraSynapses,ABD.ABDcContraSynapses];
ABD.motorRestSynapseNumbers = [ABD.ABDrRestSynapses,ABD.ABDcRestSynapses];



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
ABD.ABDrDMVL = ABD.Origin(ABD.ABDrLocs,1)+ABD.Origin(ABD.ABDrLocs,3);
ABD.ABDrDMVLindex =  (ABD.ABDrDMVL - min(ABD.ABDrDMVL)) ./ (max(ABD.ABDrDMVL)-min(ABD.ABDrDMVL));
ABD.ABDcDMVL = ABD.Origin(ABD.ABDcLocs,1)+ABD.Origin(ABD.ABDcLocs,3);
ABD.ABDcDMVLindex =  (ABD.ABDcDMVL - min(ABD.ABDcDMVL)) ./ (max(ABD.ABDcDMVL)-min(ABD.ABDcDMVL));


%%


[ABD.ABDrSaccadicPathLength,ABD.ABDrSaccadicGradient] = getABDgradient('ABDr',unique(vertcat(ABDr.Saccadic)),false,true);
[ABD.ABDrVestibularPathLength,ABD.ABDrVestibularGradient] =  getABDgradient('ABDr',unique(vertcat(ABDr.Vestibular)),false,true);
[ABD.ABDrIntegratorPathLength,ABD.ABDrIntegratorGradient] =  getABDgradient('ABDr',unique(vertcat(ABDr.Integrator)),false,true);
[ABD.ABDrContraPathLength,ABD.ABDrContraGradient] =  getABDgradient('ABDr',unique(vertcat(ABDr.Contra)),false,true);
[ABD.ABDrRestPathLength,ABD.ABDrRestGradient] =  getABDgradient('ABDr',unique(vertcat(ABDr.EverythingElse)),false,true);
[ABD.ABDrPutSaccPathLength,ABD.ABDrPutSaccGradient] =  getABDgradient('ABDr',ABDPutativeSaccadic.cellIDs,false,true);
[ABD.ABDrLatIntPathLength,ABD.ABDrLatIntGradient] =  getABDgradient('ABDr',lateralVSaccadic,false,true);
[ABD.ABDrABD_ABDiPathLength,ABD.ABDrABD_ABDiGradient] =  getABDgradient('ABDr',ABD_ABDi_PutativeSaccadic.cellIDs,false,true);


[ABD.ABDcSaccadicPathLength,ABD.ABDcSaccadicGradient] = getABDgradient('ABDc',unique(vertcat(ABDc.Saccadic)),false,true);
[ABD.ABDcVestibularPathLength,ABD.ABDcVestibularGradient] =  getABDgradient('ABDc',unique(vertcat(ABDc.Vestibular)),false,true);
[ABD.ABDcIntegratorPathLength,ABD.ABDcIntegratorGradient] =  getABDgradient('ABDc',unique(vertcat(ABDc.Integrator)),false,true);
[ABD.ABDcContraPathLength,ABD.ABDcContraGradient] =  getABDgradient('ABDc',unique(vertcat(ABDc.Contra)),false,true);
[ABD.ABDcRestPathLength,ABD.ABDcRestGradient] =  getABDgradient('ABDc',unique(vertcat(ABDc.EverythingElse)),false,true);
[ABD.ABDcPutSacccPathLength,ABD.ABDcPutSaccGradient] =  getABDgradient('ABDc',ABDPutativeSaccadic.cellIDs,false,true);
[ABD.ABDcLatIntPathLength,ABD.ABDcLatIntGradient] =  getABDgradient('ABDc',lateralVSaccadic,false,true);
[ABD.ABDcABD_ABDiPathLength,ABD.ABDcABD_ABDiGradient] =  getABDgradient('ABDc',ABD_ABDi_PutativeSaccadic.cellIDs,false,true);


% combine ABDr and ABDc
ABD.motorCellIDs  = [ABDr_CellIDs';ABDc_CellIDs'];
ABD.motorTotalSynapse = [ABD.Synapses(ABD.ABDrLocs);ABD.Synapses(ABD.ABDcLocs)];
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

ABD.motorPutSaccGradient = [horzcat(ABD.ABDrvols,ABD.ABDrPutSaccGradient);horzcat(ABD.ABDcvols,ABD.ABDcPutSaccGradient)];

ABD.motorlatIntGradient = [horzcat(ABD.ABDrvols,ABD.ABDrLatIntGradient);horzcat(ABD.ABDcvols,ABD.ABDcLatIntGradient)];
ABD.motorABD_ABDiGradient = [horzcat(ABD.ABDrvols,ABD.ABDrABD_ABDiGradient);horzcat(ABD.ABDcvols,ABD.ABDcABD_ABDiGradient)]

%%
% ABDIr

[ABD.ABDIrSaccadicPathLength,ABD.ABDIrSaccadicGradient] =  getABDgradient(ABDIr,unique(vertcat(ABDIr.Saccadic)),false,false);
[ABD.ABDIrVestibularPathLength,ABD.ABDIrVestibularGradient] =  getABDgradient(ABDIr,unique(vertcat(ABDIr.Vestibular)),false,false);
[ABD.ABDIrIntegratorPathLength,ABD.ABDIrIntegratorGradient] =  getABDgradient(ABDIr,unique(vertcat(ABDIr.Integrator)),false,false);
[ABD.ABDIrContraPathLength,ABD.ABDIrContraGradient] =  getABDgradient(ABDIr,unique(vertcat(ABDIr.Contra)),false,false);
[ABD.ABDIrRestPathLength,ABD.ABDIrRestGradient] =  getABDgradient(ABDIr,unique(vertcat(ABDIr.EverythingElse)),false,false);
[ABD.ABDIrPutSaccPathLength,ABD.ABDIrPutSaccGradient] =  getABDgradient(ABDIr,ABDiPutativeSaccadic.cellIDs,false,false);



[ABD.ABDIcSaccadicPathLength,ABD.ABDIcSaccadicGradient] =  getABDgradient(ABDIc,unique(vertcat(ABDIc.Saccadic)),false,false);
[ABD.ABDIcVestibularPathLength,ABD.ABDIcVestibularGradient] =  getABDgradient(ABDIc,unique(vertcat(ABDIc.Vestibular)),false,false);
[ABD.ABDIcIntegratorPathLength,ABD.ABDIcIntegratorGradient] =  getABDgradient(ABDIc,unique(vertcat(ABDIc.Integrator)),false,false);
[ABD.ABDIcContraPathLength,ABD.ABDIcContraGradient] =  getABDgradient(ABDIc,unique(vertcat(ABDIc.Contra)),false,false);
[ABD.ABDIcRestPathLength,ABD.ABDIcRestGradient] =  getABDgradient(ABDIc,unique(vertcat(ABDIc.EverythingElse)),false,false);
[ABD.ABDIcPutSaccPathLength,ABD.ABDIcPutSaccGradient] =  getABDgradient(ABDIc,ABDiPutativeSaccadic.cellIDs,false,false);


% combine ABDIr and ABDIc
ABD.InterCellIDs = [ABDIr_CellIDs';ABDIc_CellIDs'];
ABD.InterOrigins = []
%ABD.InterCellIDsSorted = intersect(ABD.sortOrder,find(setdiff(ABD.cellID,Allmotor)),'stable');

ABD.InterSaccadicGradient = [horzcat(ABD.ABDIrvols,ABD.ABDIrSaccadicGradient);horzcat(ABD.ABDIcvols,ABD.ABDIcSaccadicGradient)];
%ABD.InterSaccadicGradient = sortrows(ABD.InterSaccadicGradient,1);

ABD.InterVestibularGradient = [horzcat(ABD.ABDIrvols,ABD.ABDIrVestibularGradient);horzcat(ABD.ABDIcvols,ABD.ABDIcVestibularGradient)];
%ABD.InterVestibularGradient = sortrows(ABD.InterVestibularGradient,1);

ABD.InterIntegratorGradient = [horzcat(ABD.ABDIrvols,ABD.ABDIrIntegratorGradient);horzcat(ABD.ABDIcvols,ABD.ABDIcIntegratorGradient)];
%ABD.InterIntegratorGradient = sortrows(ABD.InterIntegratorGradient,1);

ABD.InterContraGradient = [horzcat(ABD.ABDIrvols,ABD.ABDIrContraGradient);horzcat(ABD.ABDIcvols,ABD.ABDIcContraGradient)];
%ABD.InterContraGradient = sortrows(ABD.InterContraGradient,1);

ABD.InterRestGradient = [horzcat(ABD.ABDIrvols,ABD.ABDIrRestGradient);horzcat(ABD.ABDIcvols,ABD.ABDIcRestGradient)];
%ABD.InterRestGradient = sortrows(ABD.InterRestGradient ,1);

ABD.InterPutSaccGradient = [horzcat(ABD.ABDIrvols,ABD.ABDIrPutSaccGradient);horzcat(ABD.ABDIcvols,ABD.ABDIcPutSaccGradient)];


%%
figure;

% r456 integrator bloc

subplot(4,4,1)
shadedErrorBar(0.1:0.1:1,nanmean(ABD.motorSaccadicGradient(:,2:end)),nanstd(ABD.motorSaccadicGradient(:,2:end))./sqrt(29),'lineprops',{'Color',ABDcolor,'Linewidth',2});
hold on;
shadedErrorBar(0.1:0.1:1,nanmean(ABD.InterSaccadicGradient(:,2:end)),nanstd(ABD.InterSaccadicGradient(:,2:end))./sqrt(21),'lineprops',{'Color',ABDicolor,'Linewidth',2});
axis square;
title('r456 Integrator axons');
legend({'M','I'},'Location','bestoutside');
xlabel('Norm. pathlength');
ylabel('Number of synapses');
offsetAxes(gca);

subplot(4,4,2)
shadedErrorBar(0.1:0.1:1,nanmean(ABD.motorlatIntGradient(:,2:end)),nanstd(ABD.motorlatIntGradient(:,2:end))./sqrt(29),'lineprops',{'Color',ABDcolor,'Linewidth',2});
hold on;
%shadedErrorBar(0.1:0.1:1,nanmean(ABD.InterSaccadicGradient(:,2:end)),nanstd(ABD.InterSaccadicGradient(:,2:end))./sqrt(21),'lineprops',{'Color',ABDicolor,'Linewidth',2});
axis square;
title('r56 Ventral Integrator axons');
legend({'M','I'},'Location','bestoutside');
xlabel('Norm. pathlength');
ylabel('Number of synapses');
offsetAxes(gca);



subplot(4,4,3)
shadedErrorBar(0.1:0.1:1,nanmean(ABD.motorIntegratorGradient(:,2:end)),nanstd(ABD.motorIntegratorGradient(:,2:end))./sqrt(29),'lineprops',{'Color',ABDcolor,'Linewidth',2});
hold on;
shadedErrorBar(0.1:0.1:1,nanmean(ABD.InterIntegratorGradient(:,2:end)),nanstd(ABD.InterIntegratorGradient(:,2:end))./sqrt(21),'lineprops',{'Color',ABDicolor,'Linewidth',2});
axis square;
title('r78 Integrator axons');
xlabel('Norm. pathlength');
ylabel('Number of synapses');
legend({'M','I'},'Location','bestoutside');
offsetAxes(gca);


% Sacccadic bloc
subplot(4,4,5)
shadedErrorBar(0.1:0.1:1,nanmean(ABD.motorPutSaccGradient(:,2:end)),nanstd(ABD.motorVestibularGradient(:,2:end))./sqrt(29),'lineprops',{'Color',ABDcolor,'Linewidth',2});
hold on;
shadedErrorBar(0.1:0.1:1,nanmean(ABD.InterPutSaccGradient(:,2:end)),nanstd(ABD.InterVestibularGradient(:,2:end))./sqrt(21),'lineprops',{'Color',ABDicolor,'Linewidth',2});
axis square;
title('Putative Saccadic axons');
xlabel('Norm. pathlength');
ylabel('Number of synapses');
legend({'M','I'},'Location','bestoutside');
offsetAxes(gca);


% Vestibular bloc


subplot(4,4,9)
shadedErrorBar(0.1:0.1:1,nanmean(ABD.motorVestibularGradient(:,2:end)),nanstd(ABD.motorVestibularGradient(:,2:end))./sqrt(29),'lineprops',{'Color',ABDcolor,'Linewidth',2});
hold on;
shadedErrorBar(0.1:0.1:1,nanmean(ABD.InterVestibularGradient(:,2:end)),nanstd(ABD.InterVestibularGradient(:,2:end))./sqrt(21),'lineprops',{'Color',ABDicolor,'Linewidth',2});
axis square;
title('r78 Integrator axons');
xlabel('Norm. pathlength');
ylabel('Number of synapses');
legend({'M','I'},'Location','bestoutside');
offsetAxes(gca);






% 
% figure;
% for i = 1:numel(ABD.motorCellIDs)
% subplot(5,6,i)
% plot(0.1:0.1:1,ABD.motorSaccadicGradient(i,2:end)./max(ABD.motorSaccadicGradient(i,2:end)));
% end


%% Arrange Medial to Lateral

saccadicColorMap = ['#ffffff'; '#f5f5f5'; '#eaeaea'; '#e0e0e0'; '#d6d6d6'; '#cccccc'; '#c2c2c2'; '#b8b8b8'; '#aeaeae'; '#a4a4a4'; '#9b9b9b'; '#919191'; '#888888'; '#7e7e7e'; '#757575'; '#6c6c6c'; '#636363'; '#5b5b5b'; '#525252'; '#4a4a4a'];
%saccadicColorMap = colorcet('L1','N',20,'reverse',1);
vestibuarColorMap = ['#ffffff'; '#fff9f4'; '#fff3e9'; '#ffedde'; '#ffe7d3'; '#ffe1c7'; '#ffdabc'; '#ffd4b1'; '#ffcea6'; '#ffc79a'; '#ffc18e'; '#ffba83'; '#ffb477'; '#ffad6a'; '#ffa65e'; '#ff9e50'; '#ff9742'; '#ff8f33'; '#ff8720'; '#ff7f00'];
integratorColorMap = ['#ffffff'; '#f5f8fb'; '#ecf1f8'; '#e2e9f4'; '#d9e2f0'; '#cfdbec'; '#c5d4e9'; '#bccde5'; '#b2c6e1'; '#a8c0dd'; '#9eb9da'; '#94b2d6'; '#8aabd2'; '#80a5ce'; '#769ecb'; '#6b98c7'; '#5f91c3'; '#538bbf'; '#4684bc'; '#377eb8'];

ABD.MotorOrigins = [vertcat(ABDr.Origin);vertcat(ABDc.Origin)];
ABD.InterOrigins = [vertcat(ABDIr.Origin);vertcat(ABDIc.Origin)];
ABD.TotalSynapses = [ABD.ABDrTotalSynapses;ABD.ABDcTotalSynapses];
% ABD.motorVols(15) = [];
% ABD.motorVols(29) = [];


figure;
% arranged Small to Big
%  [~,ABDMotorSmallBig] = sort(ABD.motorVols);
% 
% ABD.motorSaccadicGradient(15,:) = [];
% ABD.motorSaccadicGradient(29,:) = [];

% arranged ML+DV
 [~,ABDMotorMLDVSort] = sort(ABD.MotorOrigins(:,1) + ABD.MotorOrigins(:,3));
 [~,ABDInterMLDVsort] = sort(ABD.InterOrigins(:,1)  + ABD.InterOrigins(:,3));
 
 %[~,ABDMotorMLDVSort] = sort([ABD.ABDrDMVLindex;ABD.ABDcDMVLindex]);

subplot(1,6,1)
heatmap(ABD.motorSaccadicGradient(ABDMotorMLDVSort,2:end)./ABD.TotalSynapses(ABDMotorMLDVSort),'ColorScaling','scaledrows',...
    'Colormap',hex2rgb(integratorColorMap),'ColorbarVisible','off','YDisplayLabels',round(ABD.motorVols(ABDMotorMLDVSort)),'XDisplayLabels',0.1:0.1:1);
title('r4-6 Position');

% vestibular nuerons
subplot(1,6,2)
heatmap(ABD.motorVestibularGradient(ABDMotorMLDVSort,2:end)./ABD.TotalSynapses(ABDMotorMLDVSort),'ColorScaling','scaledrows',...
    'Colormap',hex2rgb(vestibuarColorMap),'ColorbarVisible','on','YDisplayLabels',round(ABD.motorVols(ABDMotorMLDVSort)),'XDisplayLabels',0.1:0.1:1);
title('Vestibular');

% 
% ABD.motorIntegratorGradient(15,:) = [];
% ABD.motorIntegratorGradient(29,:) = [];
 
subplot(1,6,3)

heatmap(ABD.motorIntegratorGradient(ABDMotorMLDVSort,2:end)./ABD.TotalSynapses(ABDMotorMLDVSort),'ColorScaling','scaledrows',...
    'Colormap',hex2rgb(integratorColorMap),'ColorbarVisible','off','YDisplayLabels',round(ABD.motorVols(ABDMotorMLDVSort)),'XDisplayLabels',0.1:0.1:1);
title('r7/8 integrator');


% Putative Saccadic Gradient

subplot(1,6,4)
heatmap(ABD.motorPutSaccGradient(ABDMotorMLDVSort,2:end)./ABD.TotalSynapses(ABDMotorMLDVSort),'ColorScaling','scaledrows',...
    'Colormap',hex2rgb(saccadicColorMap),'ColorbarVisible','on','YDisplayLabels',round(ABD.motorVols(ABDMotorMLDVSort)),'XDisplayLabels',0.1:0.1:1);
title('r2/3 Saccadic');


subplot(1,6,5)
heatmap(ABD.motorlatIntGradient(ABDMotorMLDVSort,2:end)./ABD.TotalSynapses(ABDMotorMLDVSort),'ColorScaling','scaledrows',...
    'Colormap',hex2rgb(integratorColorMap),'ColorbarVisible','off','YDisplayLabels',round(ABD.motorVols(ABDMotorMLDVSort)),'XDisplayLabels',0.1:0.1:1);
title('r56 integrator');

integratorSums = ABD.motorSaccadicGradient+ABD.motorIntegratorGradient+ABD.motorlatIntGradient;


subplot(1,6,6)
heatmap(integratorSums(ABDMotorMLDVSort,2:end)./ABD.TotalSynapses(ABDMotorMLDVSort),'ColorScaling','scaledrows',...
    'Colormap',hex2rgb(integratorColorMap),'ColorbarVisible','on','YDisplayLabels',round(ABD.motorVols(ABDMotorMLDVSort)),'XDisplayLabels',0.1:0.1:1);
title('all integrator');


figure;
subplot(4,4,1)
shadedErrorBar(0.1:0.1:1,nanmean(ABD.motorVestibularGradient(:,2:end)),nanstd(ABD.motorVestibularGradient(:,2:end))./sqrt(29),'lineprops',{'Color',[255,127,0]./255,'Linewidth',2});
hold on;
shadedErrorBar(0.1:0.1:1,nanmean(ABD.motorPutSaccGradient(:,2:end)),nanstd(ABD.motorVestibularGradient(:,2:end))./sqrt(29),'lineprops',{'Color',[74,74,74]./255,'Linewidth',2});
shadedErrorBar(0.1:0.1:1,nanmean(integratorSums(:,2:end)),nanstd(integratorSums(:,2:end))./sqrt(29),'lineprops',{'Color',[55,126,184]./255,'Linewidth',2});
axis square;
box off;
offsetAxes(gca);
xlabel('Norm. pathlength');
ylabel('Average number of synapses');

% plot by MLDV axis

figure;
subplot(4,4,1)
VtoSacratio = (sum(ABD.motorVestibularGradient(ABDMotorMLDVSort,2:end),2)./ABD.TotalSynapses(ABDMotorMLDVSort))...
./(sum(ABD.motorPutSaccGradient(ABDMotorMLDVSort,2:end),2)./ABD.TotalSynapses(ABDMotorMLDVSort));
%plot(sum(ABD.motorVestibularGradient(ABDMotorMLDVSort,2:end),2)./sum(ABD.motorPutSaccGradient(ABDMotorMLDVSort,2:end),2),'o','MarkerFaceColor',ABDcolor,'MarkerEdgeColor',ABDcolor);
plot(VtoSacratio,'o','MarkerFaceColor',ABDcolor,'MarkerEdgeColor','k');
VtoSacratio = VtoSacratio(~isinf(VtoSacratio));
VtoSacratio = VtoSacratio(~isnan(VtoSacratio));

hold on;
f = showfit(ezfit(1:size(VtoSacratio,1),VtoSacratio,'affine'),'fitcolor',ABDcolor,...
    'dispeqboxmode','off','dispfitlegend','on','corrcoefmode','r2','fitlineStyle',':');
text(size(VtoSacratio,1),max(VtoSacratio),sprintf('r=%0.2f',f.r));
view([90,90]);
xlabel('VL <--> DM');
ylabel('v/r23');
box off;
daspect([1,2,1]);
offsetAxes(gca);


movingAverageOfRatio = movmean(sum(ABD.motorVestibularGradient(ABDMotorMLDVSort,2:end),2)./sum(ABD.motorPutSaccGradient(ABDMotorMLDVSort,2:end),2),2,1,'omitnan');
movingAverageOfVolumes = movmean(round(ABD.motorVols(ABDMotorMLDVSort)),2,1);

subplot(4,4,2)
plot(movingAverageOfRatio,'o','MarkerFaceColor',ABDcolor,'MarkerEdgeColor',ABDcolor);
hold  on;
view([90,90]);
xlabel('VL <--> DM');
ylabel('moving average');
box off;
%daspect([1,1,1]);
offsetAxes(gca);

subplot(4,4,3)
plot(movingAverageOfVolumes,'o','MarkerFaceColor',ABDcolor,'MarkerEdgeColor',ABDcolor);
view([90,90]);
xlabel('VL <--> DM');
ylabel('moving average of soma vol (um^3)');
box off;
daspect([1,4,1]);
offsetAxes(gca);

subplot(4,4,5)
cmap = colorcet('R3','N',29);
scatter(ABD.motorVols(ABDMotorMLDVSort),sum(ABD.motorVestibularGradient(ABDMotorMLDVSort,2:end),2)./sum(ABD.motorPutSaccGradient(ABDMotorMLDVSort,2:end),2),40,...
    cmap,'filled','MarkerEdgeColor','k');
a = ABD.motorVols(ABDMotorMLDVSort);
b = sum(ABD.motorVestibularGradient(ABDMotorMLDVSort,2:end),2)./sum(ABD.motorPutSaccGradient(ABDMotorMLDVSort,2:end),2);
removeIndex = [find(isnan(b));find(isinf(b))];
c = a;
c(removeIndex) = [];
d = b;
d(removeIndex) = [];
f = showfit(ezfit(c,d,'c*x+d'),'fitcolor',ABDcolor,'dispeqboxmode','off','dispfitlegend','on','corrcoefmode','r2');
legend off;
text(max(c),max(d),sprintf('r=%0.2f',f.r));
box off;
axis square;
xlabel('Soma volume \mu^3');
ylabel('V/r23');
%offsetAxes(gca);


% figure; % plot by soma gradient
%  [~,ABDMotorSomaSort] = sort(ABD.motorVols);
%  subplot(4,4,1)
% VtoSacratio = (sum(ABD.motorVestibularGradient(ABDMotorSomaSort,2:end),2)./ABD.TotalSynapses(ABDMotorSomaSort))...
% ./(sum(ABD.motorPutSaccGradient(ABDMotorSomaSort,2:end),2)./ABD.TotalSynapses(ABDMotorSomaSort));
% %plot(sum(ABD.motorVestibularGradient(ABDMotorSomaSort,2:end),2)./sum(ABD.motorPutSaccGradient(ABDMotorSomaSort,2:end),2),'o','MarkerFaceColor',ABDcolor,'MarkerEdgeColor',ABDcolor);
% plot(VtoSacratio,'o','MarkerFaceColor',ABDcolor,'MarkerEdgeColor',ABDcolor);
% VtoSacratio = VtoSacratio(~isinf(VtoSacratio));
% VtoSacratio = VtoSacratio(~isnan(VtoSacratio));
% 
% hold on;
% f = showfit(ezfit(1:size(VtoSacratio,1),VtoSacratio,'affine'),'fitcolor',ABDcolor,'dispeqboxmode','off','dispfitlegend','on','corrcoefmode','r2')
% text(size(VtoSacratio,1),max(VtoSacratio),sprintf('r=%0.2f',f.r));
% view([90,90]);
% xlabel('VL <--> DM');
% ylabel('v/r23');
% box off;
% %daspect([1,1,1]);
% offsetAxes(gca);
% 
% 
% movingAverageOfRatio = movmean(sum(ABD.motorVestibularGradient(ABDMotorSomaSort,2:end),2)./sum(ABD.motorPutSaccGradient(ABDMotorSomaSort,2:end),2),2,1,'omitnan');
% movingAverageOfVolumes = movmean(round(ABD.motorVols(ABDMotorSomaSort)),2,1);
% 
% subplot(4,4,2)
% plot(movingAverageOfRatio,'o','MarkerFaceColor',ABDcolor,'MarkerEdgeColor',ABDcolor);
% hold  on;
% view([90,90]);
% xlabel('VL <--> DM');
% ylabel('moving average');
% box off;
% %daspect([1,1,1]);
% offsetAxes(gca);
% 
% subplot(4,4,3)
% plot(movingAverageOfVolumes,'o','MarkerFaceColor',ABDcolor,'MarkerEdgeColor',ABDcolor);
% view([90,90]);
% xlabel('VL <--> DM');
% ylabel('moving average of soma vol (um^3)');
% box off;
% daspect([1,4,1]);
% offsetAxes(gca);
% 
% subplot(4,4,5)
% cmap = colorcet('R3','N',29);
% scatter(ABD.motorVols(ABDMotorSomaSort),sum(ABD.motorVestibularGradient(ABDMotorSomaSort,2:end),2)./sum(ABD.motorPutSaccGradient(ABDMotorSomaSort,2:end),2),40,...
%     cmap,'filled','MarkerEdgeColor','k');
% a = ABD.motorVols(ABDMotorSomaSort);
% b = sum(ABD.motorVestibularGradient(ABDMotorSomaSort,2:end),2)./sum(ABD.motorPutSaccGradient(ABDMotorSomaSort,2:end),2);
% removeIndex = [find(isnan(b));find(isinf(b))];
% c = a;
% c(removeIndex) = [];
% d = b;
% d(removeIndex) = [];
% f = showfit(ezfit(c,d,'c*x+d'),'fitcolor',ABDcolor,'dispeqboxmode','off','dispfitlegend','on','corrcoefmode','r2');
% legend off;
% text(max(c),max(d),sprintf('r=%0.2f',f.r));
% box off;
% axis square;
% xlabel('Soma volume \mu^3');
% ylabel('V/r23');

 

%% plot for figure;

% % vestibular nuerons
% subplot(1,5,1)
% heatmap(ABD.motorVestibularGradient(ABDMotorMLDVSort,2:end),...
%     'Colormap',hex2rgb(vestibuarColorMap),'ColorbarVisible','off','YDisplayLabels',round(ABD.motorVols(ABDMotorMLDVSort)),'XDisplayLabels',0.1:0.1:1);
% title('Vestibular');
% 
% % Putative Saccadic Gradient
% 
% subplot(1,5,2)
% heatmap(ABD.motorPutSaccGradient(ABDMotorMLDVSort,2:end),...
%     'Colormap',hex2rgb(saccadicColorMap),'ColorbarVisible','off','YDisplayLabels',round(ABD.motorVols(ABDMotorMLDVSort)),'XDisplayLabels',0.1:0.1:1);
% title('r2/3 Saccadic');


% vestibular nuerons
subplot(1,6,[1])

heatmap(ABD.motorVestibularGradient(ABDMotorMLDVSort,2:end)./ABD.motorTotalSynapse(ABDMotorMLDVSort),...
    'Colormap',hex2rgb(vestibuarColorMap),'ColorbarVisible','on','XDisplayLabels',0.1:0.1:1,'ColorScaling','scaledrows','MissingDataColor','w',...
'colorlimits',[0,1]);
title('Vestibular');




% Putative Saccadic Gradient

subplot(1,6,[3])
heatmap(ABD.motorPutSaccGradient(ABDMotorMLDVSort,2:end)./ABD.motorTotalSynapse(ABDMotorMLDVSort),...
    'Colormap',hex2rgb(saccadicColorMap),'ColorbarVisible','on','XDisplayLabels',0.1:0.1:1,'ColorScaling','scaledrows','MissingDataColor','w');
title('r2/3 Saccadic');



subplot(1,6,[5])
heatmap(integratorSums(ABDMotorMLDVSort,2:end)./ABD.motorTotalSynapse(ABDMotorMLDVSort),...
    'Colormap',hex2rgb(integratorColorMap),'ColorbarVisible','on','XDisplayLabels',0.1:0.1:1,'ColorScaling','scaledrows','MissingDataColor','w');
title('all integrator');

figure;

subplot(1,6,2)
bar(0.5:1:28.5,sum(ABD.motorVestibularGradient(ABDMotorMLDVSort,2:end)./ABD.motorTotalSynapse(ABDMotorMLDVSort),2),...
    'FaceColor',[255,127,0]./255,'EdgeColor','none');
view([90,90]);
set(gca,'XLim',[0,29],'YLim',[0,0.05],'XColor','none');
ylabel('Fractional input');
box off;

subplot(1,6,4)
bar(sum(ABD.motorPutSaccGradient(ABDMotorMLDVSort,2:end)./ABD.motorTotalSynapse(ABDMotorMLDVSort),2),...
    'FaceColor',[74,74,74]./255,'EdgeColor','none');
view([90,90]);
set(gca,'XLim',[0.5,29.5],'YLim',[0,0.05],'XColor','none');
ylabel('Fractional input');
box off;

subplot(1,6,6)
bar(sum(integratorSums(ABDMotorMLDVSort,2:end)./ABD.motorTotalSynapse(ABDMotorMLDVSort),2),...
    'FaceColor',[55,126,184]./255,'EdgeColor','none');
view([90,90]);
set(gca,'XLim',[0.5,29.5],'XColor','none');
ylabel('Fractional input');
box off;


%% cumulative plots

figure;
subplot(4,4,1)
histogram([cell2mat(ABD.ABDrPutSaccPathLength),cell2mat(ABD.ABDcPutSacccPathLength)],20,'Normalization','cdf',...
    'EdgeColor',[74,74,74]./255,'LineWidth',2,'DisplayStyle','stairs');
hold on
histogram([cell2mat(ABD.ABDrVestibularPathLength),cell2mat(ABD.ABDcVestibularPathLength)],20,'Normalization','cdf',...
    'EdgeColor','r','LineWidth',2,'DisplayStyle','stairs');

IntegratorPathLengthsABD = [cell2mat(ABD.ABDrSaccadicPathLength),cell2mat(ABD.ABDcSaccadicPathLength),...
                        cell2mat(ABD.ABDrIntegratorPathLength),cell2mat(ABD.ABDcIntegratorPathLength),...
                        cell2mat(ABD.ABDrLatIntPathLength), cell2mat(ABD.ABDcLatIntPathLength)];

histogram(IntegratorPathLengthsABD,20,'Normalization','cdf',...
    'EdgeColor','b','LineWidth',2,'DisplayStyle','stairs');

axis square;
box off;
xlabel('Norm. pathlength');
ylabel('Cumulative count');
offsetAxes(gca);

subplot(4,4,2)
histogram([cell2mat(ABD.ABDIrPutSaccPathLength),cell2mat(ABD.ABDIcPutSaccPathLength)],20,'Normalization','cdf',...
    'EdgeColor',[74,74,74]./255,'LineWidth',2,'DisplayStyle','stairs');
hold on
histogram([cell2mat(ABD.ABDIrVestibularPathLength),cell2mat(ABD.ABDIcVestibularPathLength)],20,'Normalization','cdf',...
    'EdgeColor','r','LineWidth',2,'DisplayStyle','stairs');

IntegratorPathLengthsABDI = [cell2mat(ABD.ABDIrSaccadicPathLength),cell2mat(ABD.ABDIcSaccadicPathLength),...
                        cell2mat(ABD.ABDIrIntegratorPathLength),cell2mat(ABD.ABDIcIntegratorPathLength),...
                       ];

histogram(IntegratorPathLengthsABDI,20,'Normalization','cdf',...
    'EdgeColor','b','LineWidth',2,'DisplayStyle','stairs');
axis square;
box off;

legend({'Saccadic','Vestibular','Integrator'});
offsetAxes(gca);

%% sorted by soma size

[~,ABDMotorMLDVSort] = sort(ABD.motorVols);
% vestibular nuerons
subplot(1,5,1)
heatmap(ABD.motorVestibularGradient(ABDMotorMLDVSort,2:end),...
    'Colormap',hex2rgb(vestibuarColorMap),'ColorbarVisible','off','YDisplayLabels',round(ABD.motorVols(ABDMotorMLDVSort)),'XDisplayLabels',0.1:0.1:1);
title('Vestibular');

% Putative Saccadic Gradient

subplot(1,5,2)
heatmap(ABD.motorPutSaccGradient(ABDMotorMLDVSort,2:end),...
    'Colormap',hex2rgb(saccadicColorMap),'ColorbarVisible','off','YDisplayLabels',round(ABD.motorVols(ABDMotorMLDVSort)),'XDisplayLabels',0.1:0.1:1);
title('r2/3 Saccadic');


% vestibular nuerons
subplot(1,5,4)
heatmap(ABD.motorVestibularGradient(ABDMotorMLDVSort,2:end)./ABD.motorTotalSynapse(ABDMotorMLDVSort),...
    'Colormap',hex2rgb(vestibuarColorMap),'ColorbarVisible','off','YDisplayLabels',round(ABD.motorVols(ABDMotorMLDVSort)),'XDisplayLabels',0.1:0.1:1);
title('Vestibular');

% Putative Saccadic Gradient

subplot(1,5,5)
heatmap(ABD.motorPutSaccGradient(ABDMotorMLDVSort,2:end)./ABD.motorTotalSynapse(ABDMotorMLDVSort),...
    'Colormap',hex2rgb(saccadicColorMap),'ColorbarVisible','off','YDisplayLabels',round(ABD.motorVols(ABDMotorMLDVSort)),'XDisplayLabels',0.1:0.1:1);
title('r2/3 Saccadic');



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





