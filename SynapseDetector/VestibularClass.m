% Classes by anatomy
%clear;
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

%% Load Vestibular Neurons

 vestibularCellIds = [76631,76700,76656,76655,76660,76679,77393,77395,77375,...
     77815,77803,77807,80591,80579,77255,77697,77364,80223,81126,80760,81117,...
     80672,80622,80572,80569,81122,78426,81328,81575,81582,81870,82161,82165,81553,82169];
 
 MVNs = [76664 76784 76783 77394 77326 77772 77766 77756 77753 77767 77775 ...
     77780 80883 80247 80669 80698 80762 80748 80696 80761 80749 81070 81044 ...
     80991 80967 81008 80883 81045 80986 81023 78647 78619 80247 81141];
 
VestibularAxons = [vestibularCellIds,MVNs];
%VestibularAxons = [MVNs];

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

motorDistribution = vertcat(Vest.MotorDist);
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
histogram(normalizedMotorCounts,-1:0.1:1,'FaceColor','k');
hold on;
line([0.5,0.5],[0,40],'color','k','LineStyle',':');
line([-0.5,-0.5],[0,40],'color','k','LineStyle',':');
xlabel('(A-Ai)/A+Ai)');
axis square;
box off;

locs = [vertcat(ABDVestibularpop.Origin);vertcat(ABDiVestibularpop.Origin)];
groups = [ones(size(vertcat(ABDVestibularpop.Origin),1),1);2*ones(size(vertcat(ABDiVestibularpop.Origin),1),1)];

figure;
subplot(4,4,[3,4])
scatterhist(locs(:,1),locs(:,2),'Group',groups,'color',[lightRed;lightBlue]);
set(gca,'YDir','reverse');

%% Breakdown by individual populations

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
xlabel('(A-Ai)/A+Ai)');
axis square;
box off;

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
