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
    df = readtable('/Users/ashwin/Documents/SynapseDetector/04152019.csv');
else
    addpath(genpath('/usr/people/ashwinv/seungmount/research/Ashwin/Scripts'));
    df = readtable('/usr/people/ashwinv/seungmount/research/Ashwin/SynapseDetector/11252018.csv');
end

%% Load Abd Neurons

load ABDr.mat
load ABDc.mat
load ABDIr.mat
load ABDIc.mat
load ABDVols.mat

ABDr_CellIDs = [82140,82145,82143,77648,82146,77302,77301, 82194,77705,77710,77305,77672,77300,77709];
ABDc_CellIDs = [82213,77654,77646,82212,77295,81172,77657,82197,77682,77154, 77652,77658,77628,77292,77688,82195,77296];
ABDIr_CellIDs = [78553,77150,77886,78547,77631,77158,77665,78552,77668,77618,77634];
ABDIc_CellIDs = [78574,79051,77148,77625,77144,77692,77641,79066,77640,77643];

ABDvols = sortrows(ABDvols,2);
ABDvols(:,2) = ABDvols(:,2)./1e9;

ABDr_vols = ABDvols(find(ismember(ABDvols(:,1),ABDr_CellIDs)),2);
ABDc_vols = ABDvols(find(ismember(ABDvols(:,1),ABDc_CellIDs)),2);
ABDIr_vols = ABDvols(find(ismember(ABDvols(:,1),ABDIr_CellIDs)),2);
ABDIc_vols = ABDvols(find(ismember(ABDvols(:,1),ABDIc_CellIDs)),2);



%%
ContraAxons = unique([vertcat(ABDr.Contra);vertcat(ABDc.Contra);vertcat(ABDIr.Contra);vertcat(ABDIc.Contra)]);

for i = 164:numel(ContraAxons)
    Contra(i) = InputsByClass(ContraAxons(i),df,2);
end


%%

motorDistribution = vertcat(Contra.MotorDist);
noSynapsesOnMotor = find(sum(motorDistribution(:,2:5),2)==0);
motorDistribution(noSynapsesOnMotor,:) = [];
motorDistributionABD = sum(motorDistribution(:,2:3),2);
motorDistributionABDi = sum(motorDistribution(:,4:5),2);

normalizedMotorCounts = (motorDistributionABD-motorDistributionABDi)./ (motorDistributionABD+motorDistributionABDi);

ix1 = 1;
ix2 = 1;
for i = 1:size(motorDistribution,1)
    if (motorDistributionABD(i)-motorDistributionABDi(i))>0
    ABDContrapop(ix1).cellIDs = motorDistribution(i,1);
    ABDContrapop(ix1).Origin = Contra(find(ABDContrapop(ix1).cellIDs==ContraAxons)).Origin;
    ABDContrapop(ix1).MotorCounts = normalizedMotorCounts(i);
    ix1 = ix1+1;
    elseif (motorDistributionABD(i)-motorDistributionABDi(i))<0
     ABDiContrapop(ix2).cellIDs = motorDistribution(i,1);
     ABDiContrapop(ix2).Origin = Contra(find(ABDiContrapop(ix2).cellIDs==ContraAxons)).Origin;
     ABDiContrapop(ix2).MotorCounts = normalizedMotorCounts(i);
     ix2 = ix2+1;
    end
end
%%
figure;
subplot(1,2,1)
[~,index]= sort([ABDContrapop.MotorCounts]);
cmapReds = colorcet('D1','N', 2*length(index));
transform_swc_AV([ABDContrapop(index).cellIDs],cmapReds(length(index):end,:),[],true,false);
colormap(gca,cmapReds(length(index):end,:))
colorbar(gca);
subplot(1,2,2)
clear index;
[~,index]= sort([ABDiContrapop.MotorCounts]);
cmapBlues = colorcet('D1','N', 2*length(index));
transform_swc_AV([ABDiContrapop(index).cellIDs],flipud(cmapBlues(1:length(index),:)),[],true,false);
colormap(gca,flipud(cmapBlues(1:length(index),:)));
colorbar(gca);

figure;
subplot(4,4,1)
histogram(normalizedMotorCounts,-1:0.1:1,'FaceColor','k');
xlabel('(A-Ai)/A+Ai)');
axis square;
box off;

locs = [vertcat(ABDContrapop.Origin);vertcat(ABDiContrapop.Origin)];
groups = [ones(size(vertcat(ABDContrapop.Origin),1),1);2*ones(size(vertcat(ABDiContrapop.Origin),1),1)];

figure;
subplot(4,4,[3,4])
scatterhist(locs(:,1),locs(:,2),'Group',groups,'color',[lightRed;lightBlue]);
set(gca,'YDir','reverse');


%% Saccadic Motor partners

for i = 1:numel([ABDContrapop.cellIDs])
  ABDContrapop(i).MotorNeuronCounts = isPostSynapseMotor(ABDContrapop(i).cellIDs,df); 
end

for i = 1:numel([ABDiContrapop.cellIDs])
  ABDiContrapop(i).MotorNeuronCounts = isPostSynapseMotor(ABDiContrapop(i).cellIDs,df); 
end

for i = 1:numel([ABDContrapop.cellIDs])
    ABDrContrapopMotorConn(i,:) = ABDContrapop(i).MotorNeuronCounts(1,1:14); % hard coarded for number of ABDr cells
    ABDcContrapopMotorConn(i,:) = ABDContrapop(i).MotorNeuronCounts(2,1:17);
    ABDIrContrapopMotorConn(i,:) = ABDContrapop(i).MotorNeuronCounts(3,1:11);
    ABDIcContrapopMotorConn(i,:) = ABDContrapop(i).MotorNeuronCounts(4,1:10);
end

for i = 1:numel([ABDiContrapop.cellIDs])
    ABDrContrapopMotorConn(i+numel([ABDContrapop.cellIDs]),:) = ABDiContrapop(i).MotorNeuronCounts(1,1:14);
    ABDcContrapopMotorConn(i+numel([ABDContrapop.cellIDs]),:) = ABDiContrapop(i).MotorNeuronCounts(2,1:17);
    ABDIrContrapopMotorConn(i+numel([ABDContrapop.cellIDs]),:) = ABDiContrapop(i).MotorNeuronCounts(3,1:11);
    ABDIcContrapopMotorConn(i+numel([ABDContrapop.cellIDs]),:) = ABDiContrapop(i).MotorNeuronCounts(4,1:10);
end


subplot(4,4,[1,5,9])
imagesc(ABDrContrapopMotorConn);
colormap(gca,colorcet('L8'));
colorbar(gca);

subplot(4,4,[2,6,10])
imagesc(ABDcContrapopMotorConn);
colormap(gca,colorcet('L8'));
colorbar(gca);

subplot(4,4,[3,7,11])
imagesc(ABDIrContrapopMotorConn);
colormap(gca,colorcet('L8'));
colorbar(gca);

subplot(4,4,[4,8,12])
imagesc(ABDIcContrapopMotorConn);
colormap(gca,colorcet('L8'));
colorbar(gca);




for i = 1:size(ABDrContrapopMotorConn,2)
    ABDrVestCounts{i} = ABDrContrapopMotorConn((ABDrContrapopMotorConn(:,i)>0),i);
end

for i = 1:size(ABDcContrapopMotorConn,2)
    ABDcVestCounts{i} = ABDcContrapopMotorConn((ABDcContrapopMotorConn(:,i)>0),i);
end

for i = 1:size(ABDIrContrapopMotorConn,2)
    ABDIrVestCounts{i} = ABDIrContrapopMotorConn((ABDIrContrapopMotorConn(:,i)>0),i);
end

for i = 1:size(ABDIcContrapopMotorConn,2)
    ABDIcVestCounts{i} = ABDIcContrapopMotorConn((ABDIcContrapopMotorConn(:,i)>0),i);
end


subplot(4,4,13)
plot(ABDr_vols,sum(ABDrContrapopMotorConn),'o','MarkerFaceColor','k','MarkerEdgeColor','none');
box off;
subplot(4,4,14)
plot(ABDc_vols,sum(ABDcContrapopMotorConn),'o','MarkerFaceColor','k','MarkerEdgeColor','none');
box off;
subplot(4,4,15)
plot(ABDIr_vols,sum(ABDIrContrapopMotorConn),'o','MarkerFaceColor','k','MarkerEdgeColor','none');
box off;
subplot(4,4,16)
plot(ABDIc_vols,sum(ABDIcContrapopMotorConn),'o','MarkerFaceColor','k','MarkerEdgeColor','none');
box off;


%% Arrange contra by rhombomere

for i = 1:numel(ContraAxons)
    if ~isempty(Contra(i).Origin)
    ContraOrigins(i,:) = Contra(i).Origin;
    end
end

for i = 1:numel(ContraAxons)
    if ~isempty(Contra(i).Rhombomere)
        ContraRhombomeres(i,:) = Contra(i).Rhombomere(5:end);
    else
        ContraRhombomeres(i,:) = NaN;
    end
end

ContraSortRhombomere = [find(ContraRhombomeres(:,1) ==1);...
                        find(ContraRhombomeres(:,2) ==1);...
                        find(ContraRhombomeres(:,3) ==1);...
                        find(ContraRhombomeres(:,4) ==1);...
                        find(ContraRhombomeres(:,5) ==1)];

[contraOriginSorted,sortOrder] = sortrows( ContraOrigins,2);


for i = 1:numel(ContraSortRhombomere)
 ContraAxonsMotorNeuronCounts{i} = isPostSynapseMotor(ContraAxons(ContraSortRhombomere(i)),df); 
end

for i = 1:numel(ContraAxonsMotorNeuronCounts)
    ABDrSortedContra(i,:) = ContraAxonsMotorNeuronCounts{i}(1,1:14);
    ABDcSortedContra(i,:) = ContraAxonsMotorNeuronCounts{i}(2,1:17);
    ABDIrSortedContra(i,:) = ContraAxonsMotorNeuronCounts{i}(3,1:11);
    ABDIcSortedContra(i,:) = ContraAxonsMotorNeuronCounts{i}(4,1:10);
end




subplot(4,4,[1,5,9])
imagesc(ABDrSortedContra);
colormap(gca,colorcet('L8'));
colorbar(gca);

subplot(4,4,[2,6,10])
imagesc(ABDcSortedContra);
colormap(gca,colorcet('L8'));
colorbar(gca);

subplot(4,4,[3,7,11])
imagesc(ABDIrSortedContra);
colormap(gca,colorcet('L8'));
colorbar(gca);

subplot(4,4,[4,8,12])
imagesc(ABDIcSortedContra);
colormap(gca,colorcet('L8'));
colorbar(gca);

subplot(4,4,13)
plot(ABDr_vols,sum(ABDrSortedContra),'o','MarkerFaceColor','k','MarkerEdgeColor','none');
box off;
subplot(4,4,14)
plot(ABDc_vols,sum(ABDcSortedContra),'o','MarkerFaceColor','k','MarkerEdgeColor','none');
box off;
subplot(4,4,15)
plot(ABDIr_vols,sum(ABDIrSortedContra),'o','MarkerFaceColor','k','MarkerEdgeColor','none');
box off;
subplot(4,4,16)
plot(ABDIc_vols,sum(ABDIcSortedContra),'o','MarkerFaceColor','k','MarkerEdgeColor','none');
box off;


    



