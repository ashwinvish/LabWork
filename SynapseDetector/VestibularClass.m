% Classes by anatomy
clear;
addpath(genpath('/Users/ashwin/Documents/'));

colors = cbrewer('qual','Dark2',10);

temp1 = colorcet('CBTL1','N',5); %  linear-tritanopic_krjcw_5-98_c46_n256
temp2 = colorcet('CBL2','N',5);  %  linear-protanopic-deuteranopic_kbw_5-98_c40_n256

leadColor = temp1(3,:); % dark red, CB safe
lagColor = temp2(3,:); % dark blue, CB safe

PartnerColors = colorcet('CBTD1','N',10);
lightRed = PartnerColors(10,:);
lightBlue = PartnerColors(1,:);


startup

if ismac
    addpath(genpath('/Users/ashwin/Documents/'));
    df = readtable('/Users/ashwin/Google Drive/Zfish/SynapseDetector/04152019.csv');
else
    addpath(genpath('/usr/people/ashwinv/seungmount/research/Ashwin/Scripts'));
    df = readtable('/usr/people/ashwinv/seungmount/research/Ashwin/SynapseDetector/11252018.csv');
end

colorSchemes
%% Load Vestibular Neurons

vestibularCellIds = [76631,76700,76656,76655,76660,76679,77393,77395,77375,...
     77815,77803,77807,80591,80579,77255,77697,77364,80223,81126,80760,81117,...
     80672,80622,80572,80569,81122,78426,81328,81575,81582,81870,82161,82165,81553,82169];
 
 MVNs = [76664 76784 76783 77394 77326 77772 77766 77756 77753 77767 77775 ...
     77780 80883 80247 80669 80698 80762 80748 80696 80761 80749 81070 81044 ...
     80991 80967 81008 80883 81045 80986 81023 78647 78619 80247 81141];
 
 TVNs = [78046,78049,78048,78047,78050,78051,78045,80345];
 
 VestNerve = [78277,78313,80345,78238,80493,78301,82283];
 contraTVN = [76257,77743,80741,81569,80738];
 
VestibularAxons = [vestibularCellIds,MVNs,contraTVN];


for i = 1:size(vestibularCellIds,2)
    temp = SynapticPartners(vestibularCellIds(i),1,df);
    a(i) = length(temp);
    temp = temp(temp<1e5);
    temp(find(temp == vestibularCellIds(i))) = [];
    b(i) = sum(ismember(temp,vestibularCellIds));
    clear temp;
end

allvest2allvest = mean(b./a);
%VestibularAxons = [MVNs];
%%
ABDr_CellIDs = [82140,82145,82143,77648,82146,77302,77301, 82194,77705,77710,77305,77672,77300,77709];
ABDc_CellIDs = [82213,77654,77646,82212,77295,81172,77657,82197,77682,77154, 77652,77658,77628,77292,77688,82195,77296];
ABDIr_CellIDs = [78553,77150,77886,78547,77631,77158,77665,78552,77668,77618,77634];
ABDIc_CellIDs = [78574,79051,77148,77625,77144,77692,77641,79066,77640,77643];

load ABDVols.mat
ABDvols = sortrows(ABDvols,2);
ABDvols(:,2) = ABDvols(:,2)./1e9;

ABDr_vols = ABDvols(find(ismember(ABDvols(:,1),ABDr_CellIDs)),2);
ABDc_vols = ABDvols(find(ismember(ABDvols(:,1),ABDc_CellIDs)),2);
ABDIr_vols = ABDvols(find(ismember(ABDvols(:,1),ABDIr_CellIDs)),2);
ABDIc_vols = ABDvols(find(ismember(ABDvols(:,1),ABDIc_CellIDs)),2);


DONMotorPattern = isMotor(vestibularCellIds',df);
MVNmotorPattern = isMotor(MVNs',df);
TONmotorPattern = isMotor(contraTVN',df);

vestSize = size(sum(DONMotorPattern(:,2:3),2),1);
MVNsize  = size(sum(MVNmotorPattern(:,2:3),2),1);
TONsize = size(sum(TONmotorPattern(:,2:3),2),1);

figure;
subplot(4,4,1)
scatter(sum(DONMotorPattern(:,2:3),2),sum(DONMotorPattern(:,4:5),2),40,Vestcolor,'filled','MarkerEdgeColor','k');
hold on,
scatter(sum(MVNmotorPattern(:,2:3),2),sum(MVNmotorPattern(:,4:5),2),40,MVNcolor,'filled','MarkerEdgeColor','k');
scatter(sum(TONmotorPattern(:,2:3),2),sum(TONmotorPattern(:,4:5),2),40,TONcolor,'filled','MarkerEdgeColor','k');

axis square;
box off;
set(gca,'XLim',[0,50],'YLim',[0,50]);
xlabel('Synapses onto M');
ylabel('Synapses onto I');
offsetAxes;



figure('Position',[100 100 350 300]);
g = gramm('x',sum(DONMotorPattern(:,2:3),2),...
    'y',sum(DONMotorPattern(:,4:5),2)); 
g.set_color_options('map',Vestcolor);%[repmat(Vestcolor,vestSize,1);repmat(MVNcolor,MVNsize,1)],'n_color',2,'n_lightness',1);
g.geom_point();
g.geom_abline();
g.stat_cornerhist('location',50,'edges',-50:5:50,'aspect',0.4);
%g(1,1).set_limit_extra( [0.05,0.05],[0.05,0.05]);
g.axe_property('XLim',[0,80],'YLim',[0,80]);
g.set_text_options('base_size',15);
g.set_names('x','Synapses on M','y','Synapses on I');
g.set_title('Vestibular')
%g.draw();

g.update('x',sum(MVNmotorPattern(:,2:3),2),'y',sum(MVNmotorPattern(:,4:5),2));
g.set_color_options('map',MVNcolor);
g.geom_point();
g.geom_abline();
g.stat_cornerhist('location',50,'edges',-50:5:50, 'aspect',0.4);
g.axe_property('XLim',[0,80],'YLim',[0,80]);
g.draw();
g.export('file_name','VestibularPlot','export_path','/Users/ashwin/Desktop','file_type','svg')

%%
 for i = 1:numel(VestibularAxons)
     if ~isExistReRoot(VestibularAxons(i))==0
     Vest(i) = InputsByClass(VestibularAxons(i),df,2);
     end
 end
 
%% Fraction reconstructed

for i= 1:numel(Vest)
    if ~isempty(Vest(i).cellID)
    [a,~] = SynapticPartners(Vest(i).cellID,1,df);
    NumberOfPrePartners(i) = numel(a);
    NumberOfPrePartnersReconstructed(i) = numel(a(a<1e5));
    clear a;
    end
end

%%

motorDistribution = isMotor(vertcat(Vest.MotorDist),df);
noSynapsesOnMotor = find(sum(motorDistribution(:,2:5),2)==0);
motorDistribution(noSynapsesOnMotor,:) = [];
motorDistributionABD = sum(motorDistribution(:,2:3),2);
motorDistributionABDi = sum(motorDistribution(:,4:5),2);

normalizedMotorCounts = (motorDistributionABD-motorDistributionABDi)./ (motorDistributionABD+motorDistributionABDi);

ix1 = 1;
ix2 = 1;
for i = 1:size(motorDistribution,1)
    if (motorDistributionABD(i)-motorDistributionABDi(i))>0
    ABDVestibularpop(ix1).cellIDs = motorDistribution(i,1);
    ABDVestibularpop(ix1).Origin = Vest(find(ABDVestibularpop(ix1).cellIDs==VestibularAxons)).Origin;
    ABDVestibularpop(ix1).MotorCounts = normalizedMotorCounts(i);
    ix1 = ix1+1;
    elseif (motorDistributionABD(i)-motorDistributionABDi(i))<0
     ABDiVestibularpop(ix2).cellIDs = motorDistribution(i,1);
     ABDiVestibularpop(ix2).Origin = Vest(find(ABDiVestibularpop(ix2).cellIDs==VestibularAxons)).Origin;
     ABDiVestibularpop(ix2).MotorCounts = normalizedMotorCounts(i);
     ix2 = ix2+1;
    end
end
%%
figure;
subplot(1,2,1)
[~,index]= sort([ABDVestibularpop.MotorCounts]);
cmapReds = colorcet('D1','N', 2*length(index));
transform_swc_AV([ABDVestibularpop(index).cellIDs],cmapReds(length(index):end,:),[],true,false);
colormap(gca,cmapReds(length(index):end,:))
colorbar(gca);
subplot(1,2,2)
clear index;
[~,index]= sort([ABDiVestibularpop.MotorCounts]);
cmapBlues = colorcet('D1','N', 2*length(index));
transform_swc_AV([ABDiVestibularpop(index).cellIDs],flipud(cmapBlues(1:length(index),:)),[],true,false);
colormap(gca,flipud(cmapBlues(1:length(index),:)));
colorbar(gca);

figure;
subplot(4,4,1)
histogram(normalizedMotorCounts(1:35),-1:0.1:1,'FaceColor',Vestcolor);
hold on;
histogram(normalizedMotorCounts(36:end-5),-1:0.1:1,'FaceColor',MVNcolor);
histogram(normalizedMotorCounts(end-5:end),-1:0.1:1,'FaceColor',TONcolor);
% line([0.5,0.5],[0,10],'color','k','LineStyle',':');
% line([-0.5,-0.5],[0,10],'color','k','LineStyle',':');
xlabel('(A-Ai)/A+Ai)');
set(gca,'YLim',[0,10]);
axis square;
box off;
offsetAxes(gca);

locs = [vertcat(ABDVestibularpop.Origin);vertcat(ABDiVestibularpop.Origin)];
groups = [ones(size(vertcat(ABDVestibularpop.Origin),1),1);2*ones(size(vertcat(ABDiVestibularpop.Origin),1),1)];

figure;
subplot(4,4,[3,4])
scatterhist(locs(:,1),locs(:,2),'Group',groups,'color',[lightRed;lightBlue]);
set(gca,'YDir','reverse');

%% Breakdown by individual populations
figure;
transform_swc_AV(vestibularCellIds,Vestcolor,[],true,false);
figure;
transform_swc_AV(MVNs,MVNcolor,[],true,false);

numberOfVestCells = numel(vestibularCellIds);

figure;
subplot(4,4,1)
histogram(normalizedMotorCounts(1:numberOfVestCells),-1:0.1:1,'FaceColor',Vestcolor);
hold on;
histogram(normalizedMotorCounts(numberOfVestCells+1:end),-1:0.1:1,'FaceColor',MVNcolor);

line([0.5,0.5],[0,10],'color','k','LineStyle',':');
line([-0.5,-0.5],[0,10],'color','k','LineStyle',':');
xlabel('OSI');
axis square;
box off;
offsetAxes

%% Vestibular Motor partners

for i = 1:numel([ABDVestibularpop.cellIDs])
  ABDVestibularpop(i).MotorNeuronCounts = isPostSynapseMotor(ABDVestibularpop(i).cellIDs,df); 
end

for i = 1:numel([ABDiVestibularpop.cellIDs])
  ABDiVestibularpop(i).MotorNeuronCounts = isPostSynapseMotor(ABDiVestibularpop(i).cellIDs,df); 
end

for i = 1:numel([ABDVestibularpop.cellIDs])
    ABDrVestibularpopMotorConn(i,:) = ABDVestibularpop(i).MotorNeuronCounts(1,1:14); % hard coarded for number of ABDr cells
    ABDcVestibularpopMotorConn(i,:) = ABDVestibularpop(i).MotorNeuronCounts(2,1:17);
    ABDIrVestibularpopMotorConn(i,:) = ABDVestibularpop(i).MotorNeuronCounts(3,1:11);
    ABDIcVestibularpopMotorConn(i,:) = ABDVestibularpop(i).MotorNeuronCounts(4,1:10);
end
ABDvestPop = [ABDVestibularpop.cellIDs];

for i = 1:numel([ABDiVestibularpop.cellIDs])
    ABDrVestibularpopMotorConn(i+numel([ABDVestibularpop.cellIDs]),:) = ABDiVestibularpop(i).MotorNeuronCounts(1,1:14);
    ABDcVestibularpopMotorConn(i+numel([ABDVestibularpop.cellIDs]),:) = ABDiVestibularpop(i).MotorNeuronCounts(2,1:17);
    ABDIrVestibularpopMotorConn(i+numel([ABDVestibularpop.cellIDs]),:) = ABDiVestibularpop(i).MotorNeuronCounts(3,1:11);
    ABDIcVestibularpopMotorConn(i+numel([ABDVestibularpop.cellIDs]),:) = ABDiVestibularpop(i).MotorNeuronCounts(4,1:10);
end

ABDiVestPop = [ABDiVestibularpop.cellIDs];

subplot(4,4,[1,5,9])
imagesc(ABDrVestibularpopMotorConn);
colormap(gca,colorcet('L8'));
colorbar(gca);

subplot(4,4,[2,6,10])
imagesc(ABDcVestibularpopMotorConn);
colormap(gca,colorcet('L8'));
colorbar(gca);

subplot(4,4,[3,7,11])
imagesc(ABDIrVestibularpopMotorConn);
colormap(gca,colorcet('L8'));
colorbar(gca);

subplot(4,4,[4,8,12])
imagesc(ABDIcVestibularpopMotorConn);
colormap(gca,colorcet('L8'));
colorbar(gca);




for i = 1:size(ABDrVestibularpopMotorConn,2)
    ABDrVestCounts{i} = ABDrVestibularpopMotorConn((ABDrVestibularpopMotorConn(:,i)>0),i);
end

for i = 1:size(ABDcVestibularpopMotorConn,2)
    ABDcVestCounts{i} = ABDcVestibularpopMotorConn((ABDcVestibularpopMotorConn(:,i)>0),i);
end

for i = 1:size(ABDIrVestibularpopMotorConn,2)
    ABDIrVestCounts{i} = ABDIrVestibularpopMotorConn((ABDIrVestibularpopMotorConn(:,i)>0),i);
end

for i = 1:size(ABDIcVestibularpopMotorConn,2)
    ABDIcVestCounts{i} = ABDIcVestibularpopMotorConn((ABDIcVestibularpopMotorConn(:,i)>0),i);
end



subplot(4,4,13)
plot(ABDr_vols,sum(ABDrVestibularpopMotorConn),'o','MarkerFaceColor','k','MarkerEdgeColor','none');
box off;
subplot(4,4,14)
plot(ABDc_vols,sum(ABDcVestibularpopMotorConn),'o','MarkerFaceColor','k','MarkerEdgeColor','none');
box off;
subplot(4,4,15)
plot(ABDIr_vols,sum(ABDIrVestibularpopMotorConn),'o','MarkerFaceColor','k','MarkerEdgeColor','none');
box off;
subplot(4,4,16)
plot(ABDIc_vols,sum(ABDIcVestibularpopMotorConn),'o','MarkerFaceColor','k','MarkerEdgeColor','none');
box off;

%% Where on the ABD are the synapses

load ABDr.mat
load ABDc.mat
load ABDIr.mat
load ABDIc.mat

[ABDVestPathLengthABDr,ABDrVestibularGradient] = getABDgradient(ABDr,vestibularCellIds,true,false);
[ABDVestPathLengthABDc,ABDcVestibularGradient] = getABDgradient(ABDc,vestibularCellIds,true,false);

[ABDVestPathLengthABDIr,ABDIrVestibularGradient] = getABDgradient(ABDIr,vestibularCellIds,true,false);
[ABDVestPathLengthABDIc,ABDIcVestibularGradient] = getABDgradient(ABDIc,vestibularCellIds,true,false);


[ABDMVNPathLengthABDr,ABDrMVNGradient] = getABDgradient(ABDr,MVNs,true,false);
[ABDMVNPathLengthABDc,ABDcMVNGradient] = getABDgradient(ABDc,MVNs,true,false);

[ABDMVNPathLengthABDIr,ABDIrMVNGradient] = getABDgradient(ABDIr,MVNs,true,false);
[ABDMVNPathLengthABDIc,ABDIcMVNGradient] = getABDgradient(ABDIc,MVNs,true,false);


[ABDTONPathLengthABDr,ABDrTONGradient] = getABDgradient(ABDr,contraTVN,true,false);
[ABDTONPathLengthABDc,ABDcTONGradient] = getABDgradient(ABDc,contraTVN,true,false);

[ABDTONPathLengthABDIr,ABDIrTONGradient] = getABDgradient(ABDIr,contraTVN,true,false);
[ABDTONPathLengthABDIc,ABDIcTONGradient] = getABDgradient(ABDIc,contraTVN,true,false);


figure;

subplot(4,4,1)
histogram([cell2mat(ABDVestPathLengthABDr),cell2mat(ABDVestPathLengthABDc)],'EdgeColor','r','DisplayStyle','stairs','LineWidth',2);
hold on;
histogram([cell2mat(ABDVestPathLengthABDIr),cell2mat(ABDVestPathLengthABDIc)],'EdgeColor','b','DisplayStyle','stairs','LineWidth',2);
axis square;
box off;
xlabel('Norm. pathlength');
ylabel('counts');
title('DO');

subplot(4,4,2)
histogram([cell2mat(ABDMVNPathLengthABDr),cell2mat(ABDMVNPathLengthABDc)],'EdgeColor','r','DisplayStyle','stairs','LineWidth',2);
hold on;
histogram([cell2mat(ABDMVNPathLengthABDIr),cell2mat(ABDMVNPathLengthABDIc)],'EdgeColor','b','DisplayStyle','stairs','LineWidth',2);
axis square;
box off;
xlabel('Norm. pathlength');
ylabel('counts');
title('MVN');



%%



% for i = 1:numel(ABDr_CellIDs)
%     temp1 = ismember(ABDr(i).Inputs,vestibularCellIds);
%     ABDVestPathLengthABDr{i} = ABDr(i).PathLength(temp1)'./max(Pvec_tree(ABDr(i).Tree{1}));
%     clear temp1;
%     
%     temp1 = ismember(ABDr(i).Inputs,MVNs);
%     ABDMVNPathLengthABDr{i} = ABDr(i).PathLength(temp1)'./max(Pvec_tree(ABDr(i).Tree{1}));
%     clear temp1;
% end
% 
% for i = 1:numel(ABDr_CellIDs)
%    ABDrVestibularGradient(i,:) =  histcounts(ABDVestPathLengthABDr{i},0:0.1:1);
%    ABDrVestibularGradient(i,:) = ABDrVestibularGradient(i,:)/max(ABDrVestibularGradient(i,:));
%    
%    ABDrMVNGradient(i,:) =  histcounts(ABDMVNPathLengthABDr{i},0:0.1:1);
%    ABDrMVNGradient(i,:) = ABDrMVNGradient(i,:)/max(ABDrMVNGradient(i,:));
% end
% 
% 
% for i = 1:numel(ABDc_CellIDs)
%     if ~isempty(ABDc(i).Tree)
%         temp1 = ismember(ABDc(i).Inputs,vestibularCellIds);
%         ABDVestPathLengthABDc{i} = ABDc(i).PathLength(temp1)'./max(Pvec_tree(ABDc(i).Tree{1}));
%         clear temp1;
%         
%         temp1 = ismember(ABDc(i).Inputs,MVNs);
%         ABDMVNPathLengthABDc{i} = ABDc(i).PathLength(temp1)'./max(Pvec_tree(ABDc(i).Tree{1}));
%         clear temp1;
%     end
% end
% 
% for i = 1:numel(ABDc_CellIDs)
%    ABDcVestibularGradient(i,:) =  histcounts(ABDVestPathLengthABDc{i},0:0.1:1);
%    ABDcVestibularGradient(i,:) = ABDcVestibularGradient(i,:)/max(ABDcVestibularGradient(i,:));
%    
%    ABDcMVNGradient(i,:) =  histcounts(ABDMVNPathLengthABDc{i},0:0.1:1);
%    ABDcMVNGradient(i,:) = ABDcMVNGradient(i,:)/max(ABDcMVNGradient(i,:));
% end
% 
% 
% for i = 1:numel(ABDIr_CellIDs)
%     temp1 = ismember(ABDIr(i).Inputs,vestibularCellIds);
%     ABDVestPathLengthABDIr{i} = ABDIr(i).PathLength(temp1)'./max(Pvec_tree(ABDIr(i).Tree{1}));
%     clear temp1;
%     
%     temp1 = ismember(ABDIr(i).Inputs,MVNs);
%     ABDMVNPathLengthABDIr{i} = ABDIr(i).PathLength(temp1)'./max(Pvec_tree(ABDIr(i).Tree{1}));
%     clear temp1;
% end
% 
% for i = 1:numel(ABDIr_CellIDs)
%    ABDIrVestibularGradient(i,:) =  histcounts(ABDVestPathLengthABDIr{i},0:0.1:1);
%    ABDIrVestibularGradient(i,:) = ABDIrVestibularGradient(i,:)/max(ABDIrVestibularGradient(i,:));
%    
%    ABDIrMVNGradient(i,:) =  histcounts(ABDMVNPathLengthABDIr{i},0:0.1:1);
%    ABDIrMVNGradient(i,:) = ABDIrMVNGradient(i,:)/max(ABDIrMVNGradient(i,:));
% end
% 
% for i = 1:numel(ABDIc_CellIDs)
%     if ~isempty(ABDIc(i).Tree)
%         temp1 = ismember(ABDIc(i).Inputs,vestibularCellIds);
%         ABDVestPathLengthABDIc{i} = ABDIc(i).PathLength(temp1)'./max(Pvec_tree(ABDIc(i).Tree{1}));
%         clear temp1;
%         
%         temp1 = ismember(ABDIc(i).Inputs,MVNs);
%         ABDMVNPathLengthABDIc{i} = ABDIc(i).PathLength(temp1)'./max(Pvec_tree(ABDIc(i).Tree{1}));
%         clear temp1;
%     end
% end
% 
% for i = 1:numel(ABDIc_CellIDs)
%    ABDIcVestibularGradient(i,:) =  histcounts(ABDVestPathLengthABDIc{i},0:0.1:1);
%    ABDIcVestibularGradient(i,:) = ABDIcVestibularGradient(i,:)/max(ABDIcVestibularGradient(i,:));
%    
%    ABDIcMVNGradient(i,:) =  histcounts(ABDMVNPathLengthABDIc{i},0:0.1:1);
%    ABDIcMVNGradient(i,:) = ABDIcMVNGradient(i,:)/max(ABDIcMVNGradient(i,:));
% end

meanVestibularABDgradient = nanmean(vertcat(ABDrVestibularGradient,ABDcVestibularGradient));
stdVestibularABDgradient = nanstd(vertcat(ABDrVestibularGradient,ABDcVestibularGradient));
numberOfABDneurons = 29;

meanVeastibularABDigradient = nanmean(vertcat(ABDIrVestibularGradient,ABDIcVestibularGradient));
stdVestibularABDigradient = nanstd(vertcat(ABDIrVestibularGradient,ABDIcVestibularGradient));
numberOfABDineurons = 21;

meanMVNABDgradient = nanmean(vertcat(ABDrMVNGradient,ABDcMVNGradient));
stdMVNABDgradient = nanstd(vertcat(ABDrMVNGradient,ABDcMVNGradient));
numberOfABDneurons = 29;

meanMVNABDigradient = nanmean(vertcat(ABDIrMVNGradient,ABDIcMVNGradient));
stdMVNABDigradient = nanstd(vertcat(ABDIrMVNGradient,ABDIcMVNGradient));
numberOfABDineurons = 21;

meanTONABDgradient = nanmean(vertcat(ABDrTONGradient,ABDcTONGradient));
stdTONABDgradient = nanstd(vertcat(ABDrTONGradient,ABDcTONGradient));
numberOfABDneurons = 29;

meanTONABDigradient = nanmean(vertcat(ABDIrTONGradient,ABDIcTONGradient));
stdTONABDigradient = nanstd(vertcat(ABDIrTONGradient,ABDIcTONGradient));
numberOfABDineurons = 21;

figure('Units','normalized','Position',[0,0,1,1]);
subplot(4,4,1)
errorbar(meanVestibularABDgradient,stdVestibularABDgradient./sqrt(numberOfABDneurons),...
    '-o','color',Vestcolor,'LineWidth',2,'MarkerFaceColor','w');
hold on;

errorbar(meanMVNABDgradient,stdMVNABDgradient./sqrt(numberOfABDneurons),...
    '-o','color',MVNcolor,'LineWidth',2,'MarkerFaceColor','w');

errorbar(meanTONABDgradient,stdTONABDgradient./sqrt(numberOfABDneurons),...
    '-o','color',TONcolor,'LineWidth',2,'MarkerFaceColor','w');

axis square;
set(gca,'XTickLabels',[0,0.5,1]);
xlabel('Norm. pathlength');
ylabel('Norm. count');
box off;
offsetAxes(gca);

subplot(4,4,2)
errorbar(meanVeastibularABDigradient,stdVestibularABDigradient./sqrt(numberOfABDineurons),...
    '-o','color',Vestcolor,'LineWidth',2,'MarkerFaceColor','w');
hold on;
errorbar(meanMVNABDigradient,stdMVNABDigradient./sqrt(numberOfABDineurons),...
    '-o','color',MVNcolor,'LineWidth',2,'MarkerFaceColor','w');
errorbar(meanTONABDigradient,stdTONABDigradient./sqrt(numberOfABDineurons),...
    '-o','color',TONcolor,'LineWidth',2,'MarkerFaceColor','w');
axis square;
set(gca,'XTickLabels',[0,0.5,1]);
xlabel('Norm. pathlength');
ylabel('Norm. count');
box off;
offsetAxes(gca);
