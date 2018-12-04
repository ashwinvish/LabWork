addpath(genpath('../'));
colors = cbrewer('qual', 'Set1', 10);
startup;

% DBX - red
% ALX - blue
% Barhl - green

if ismac
    addpath(genpath('/Users/admin/Documents/Scripts'));
    df = readtable('/Users/admin/Documents/SynapseDetector/09202018.csv');
else
    addpath(genpath('/usr/people/ashwinv/seungmount/research/Ashwin/Scripts'));
    df = readtable('/usr/people/ashwinv/seungmount/research/Ashwin/SynapseDetector/09202018.csv');
end

%load('FiringRates.mat');
load('TAU.mat');
load('STA.mat');
load('AllCells.mat');
Firing = A;

DbxCells = [76182,76183,76185,76186,76188,76189,76191,76199,76200];
OriginalCellOrderDBX = [8,9,10,11,21,12,15,2,3];
DbxTimeConstants = TAU(OriginalCellOrderDBX);

AlxCells = [76181,76201,76187,76184,76192,76197];
OriginalCellOrderALX = [5,7,6,13,16,22];
AlxTimeConstants = TAU(OriginalCellOrderALX);

BarhlCells = [76198,76190,76193,76194,76195,76196];
OriginalCellOrderBARHL = [1,14,17,18,19,20];
BarhlTimeConstants = TAU(OriginalCellOrderBARHL);

AllCellClasses = [DbxCells,AlxCells,BarhlCells];
OriginalCellOrder =  [OriginalCellOrderDBX,OriginalCellOrderALX,OriginalCellOrderBARHL];


%%
% DBX
corrAndSizeDBX = [];
for i = 1:length(DbxCells)
    for j = 1:length(DbxCells)
        if i == j
            corrAndSizeDBX = [corrAndSizeDBX];
        else
            [preSynapticInputs1,preSynapticInputs1PSD] = SynapticPartners(DbxCells(i),1,df);
            [preSynapticInputs2,preSynapticInputs2PSD] = SynapticPartners(DbxCells(j),1,df);
            [commonInputs,loc1,loc2] = intersect(preSynapticInputs1, preSynapticInputs2);
            psdSize1 = df.size(preSynapticInputs1PSD(loc1));
            psdSize2 = df.size(preSynapticInputs2PSD(loc2));
            corr = corrcoef(Firing(:,OriginalCellOrderDBX(i)), Firing(:,OriginalCellOrderDBX(j)));
            corrAndSizeDBX = [corrAndSizeDBX; commonInputs, corr(1,2).*ones(size(psdSize1,1),1), psdSize1, psdSize2];
            
            clear {'preSynapticInputs1','preSynapticInputs1PSD', 'psdSize1','loc1'};
            clear {'preSynapticInputs2', 'preSynapticInputs2PSD','psdSize2', 'loc2'};
            clear commonInputs;
        end
    end
end


% ALX
corrAndSizeALX = [];
for i = 1:length(AlxCells)
    for j = 1:length(AlxCells)
        if i == j
            corrAndSizeALX = [corrAndSizeALX];
        else
            [preSynapticInputs1,preSynapticInputs1PSD] = SynapticPartners(AlxCells(i),1,df);
            [preSynapticInputs2,preSynapticInputs2PSD] = SynapticPartners(AlxCells(j),1,df);
            [commonInputs,loc1,loc2] = intersect(preSynapticInputs1, preSynapticInputs2);
            psdSize1 = df.size(preSynapticInputs1PSD(loc1));
            psdSize2 = df.size(preSynapticInputs2PSD(loc2));
            corr = corrcoef(Firing(:,OriginalCellOrderALX(i)), Firing(:,OriginalCellOrderALX(j)));
            corrAndSizeALX = [corrAndSizeALX;commonInputs, corr(1,2).*ones(size(psdSize1,1),1), psdSize1, psdSize2];
            
            clear {'preSynapticInputs1','preSynapticInputs1PSD', 'psdSize1','loc1'};
            clear {'preSynapticInputs2', 'preSynapticInputs2PSD','psdSize2', 'loc2'};
            clear commonInputs;
        end
    end
end

%BARHL
corrAndSizeBARHL = [];
for i = 1:length(BarhlCells)
    for j = 1:length(BarhlCells)
        if i == j
            corrAndSizeBARHL = [corrAndSizeBARHL];
        else
            [preSynapticInputs1,preSynapticInputs1PSD] = SynapticPartners(BarhlCells(i),1,df);
            [preSynapticInputs2,preSynapticInputs2PSD] = SynapticPartners(BarhlCells(j),1,df);
            [commonInputs,loc1,loc2] = intersect(preSynapticInputs1, preSynapticInputs2);
            psdSize1 = df.size(preSynapticInputs1PSD(loc1));
            psdSize2 = df.size(preSynapticInputs2PSD(loc2));
            corr = corrcoef(Firing(:,OriginalCellOrderBARHL(i)), Firing(:,OriginalCellOrderBARHL(j)));
            corrAndSizeBARHL = [corrAndSizeBARHL;commonInputs, corr(1,2).*ones(size(psdSize1,1),1), psdSize1, psdSize2];
            
            clear {'preSynapticInputs1','preSynapticInputs1PSD', 'psdSize1','loc1'};
            clear {'preSynapticInputs2', 'preSynapticInputs2PSD','psdSize2', 'loc2'};
            clear commonInputs;
        end
    end
end

% All cells and all coorelations

corrAndSize = [];
for i = 1:length(AllCellClasses)
    for j = 1:length(AllCellClasses)
        if i == j
            corrAndSize = [corrAndSize];
        else
            [preSynapticInputs1,preSynapticInputs1PSD] = SynapticPartners(AllCellClasses(i),1,df);
            [preSynapticInputs2,preSynapticInputs2PSD] = SynapticPartners(AllCellClasses(j),1,df);
            [commonInputs,loc1,loc2] = intersect(preSynapticInputs1, preSynapticInputs2);
            psdSize1 = df.size(preSynapticInputs1PSD(loc1));
            psdSize2 = df.size(preSynapticInputs2PSD(loc2));
            corr = corrcoef(Firing(:,OriginalCellOrder(i)), Firing(:,OriginalCellOrder(j)));
            corrAndSize = [corrAndSize; commonInputs, corr(1,2).*ones(size(commonInputs,1),1), psdSize1, psdSize2];
            
            clear {'preSynapticInputs1','preSynapticInputs1PSD', 'psdSize1','loc1'};
            clear {'preSynapticInputs2', 'preSynapticInputs2PSD','psdSize2', 'loc2'};
            clear commonInputs;
        end
    end
end


subplot(4,4,5)
scatter( (corrAndSizeDBX(:,3)+corrAndSizeDBX(:,4)), corrAndSizeDBX(:,2), 20, 'MarkerFaceColor',colors(1,:),'MarkerEdgeColor','none');
%showfit('linear','fitcolor',colors(1,:), 'Linestyle','--');
set(gca,'YLim',[-0.3,1],'XLim', [min(corrAndSize(:,3)+corrAndSize(:,4)), max(corrAndSize(:,3)+corrAndSize(:,4))])
ylabel('STA correlation');
xlabel('PSD size (voxels)');
axis square;
offsetAxes

subplot(4,4,6)
scatter((corrAndSizeALX(:,3)+corrAndSizeALX(:,4)), corrAndSizeALX(:,2), 20, 'MarkerFaceColor',colors(2,:),'MarkerEdgeColor','none');
set(gca,'YLim',[-0.3,1])
xlabel('STA correlation');
ylabel('PSD size (voxels)');
axis square;
%showfit('linear','fitcolor',colors(2,:), 'Linestyle','--');
offsetAxes

subplot(4,4,7)
scatter((corrAndSizeBARHL(:,3)+corrAndSizeBARHL(:,4)), corrAndSizeBARHL(:,2), 20, 'MarkerFaceColor',colors(3,:),'MarkerEdgeColor','none');
set(gca,'YLim',[-0.3,1])
xlabel('STA correlation');
ylabel('PSD size (voxels)');
axis square;
%showfit('linear','fitcolor',colors(3,:), 'Linestyle','--');
offsetAxes


subplot(4,4,9)
scatter((corrAndSize(:,3)+corrAndSize(:,4)),corrAndSize(:,2) ,30, 'MarkerFaceColor',colors(4,:),'MarkerEdgeColor','none');
set(gca,'YLim',[-0.3,1])
xlabel('STA correlation');
ylabel('PSD size (voxels)');
axis square;
%showfit('linear','fitcolor',colors(4,:), 'Linestyle','--');
offsetAxes


%% Break up axon analysis by axon type

for i = 1:size(corrAndSizeDBX,1)
    


for i = 1:size(CommonInputsDBX,2)
    [a,b,c] = intersect(mlOrdered,CommonInputsDBX{i});
    commonInputsCuratedDBX{i} = a;
    commonInputsCuratedTypeDBX{i} = mlOrderedTypes(b);
end

% VestibularNeurons
for i =1:size(CommonInputsDBX,2)
   index =  find(strcmp(commonInputsCuratedTypeDBX{i},'Vestibular'));
   commonDBXVestibularInputs{i} = commonInputsCuratedDBX{i}(index);  
   clear index;
   index =  find(strcmp(commonInputsCuratedTypeDBX{i},'Bushy'));
   commonDBXBushyInputs{i} = commonInputsCuratedDBX{i}(index);  
   clear index;
   index =  find(strcmp(commonInputsCuratedTypeDBX{i},'Lateral_Dorsal'));
   commonDBXLateralDInputs{i} = commonInputsCuratedDBX{i}(index);
   clear index;
   index =  find(strcmp(commonInputsCuratedTypeDBX{i},'PutativeALX'));
   commonDBXPutALXInputs{i} = commonInputsCuratedDBX{i}(index);
end

for i = 1:length(DbxCells)
    for j = 1:length(DbxCells)
        if i==j
            sumOfCommonVestibularInputs(i,j) = NaN;
            sumOfCommonBushyInputs(i,j) = NaN;
            sumOfCommonLateralDInputs(i,j) = NaN;
            sumOfCommonPutALXInputs(i,j) = NaN;
            individualDBXVestibularInputs{sub2ind(size(FiringCorrDBX),i,j)} = [NaN, NaN];
            individualDBXBushyInputs{sub2ind(size(FiringCorrDBX),i,j)} = [NaN, NaN];
            individualDBXLateralDInputs{sub2ind(size(FiringCorrDBX),i,j)} = [NaN, NaN];
            individualDBXPutALXInputs{sub2ind(size(FiringCorrDBX),i,j)} = [NaN, NaN];
        else
            [temp1,temp2]= commonInputSize(DbxCells(i),DbxCells(j),commonDBXVestibularInputs{sub2ind(size(FiringCorrDBX),i,j)},df);
            sumOfCommonVestibularInputs(i,j) = sum([sum(temp1),sum(temp2)]);
            individualDBXVestibularInputs{sub2ind(size(FiringCorrDBX),i,j)} = [sum(temp1),sum(temp2)];
            clear {temp1,temp2};
            [temp1,temp2]= commonInputSize(DbxCells(i),DbxCells(j),commonDBXBushyInputs{sub2ind(size(FiringCorrDBX),i,j)},df);
            sumOfCommonBushyInputs(i,j) = sum([sum(temp1),sum(temp2)]);
            individualDBXBushyInputs{sub2ind(size(FiringCorrDBX),i,j)} = [sum(temp1),sum(temp2)];
            clear {temp1,temp2};
            [temp1,temp2]= commonInputSize(DbxCells(i),DbxCells(j),commonDBXLateralDInputs{sub2ind(size(FiringCorrDBX),i,j)},df);
            sumOfCommonLateralDInputs(i,j) = sum([sum(temp1),sum(temp2)]);
            individualDBXLateralDInputs{sub2ind(size(FiringCorrDBX),i,j)} = [sum(temp1),sum(temp2)];
            clear {temp1,temp2};
            [temp1,temp2]= commonInputSize(DbxCells(i),DbxCells(j),commonDBXPutALXInputs{sub2ind(size(FiringCorrDBX),i,j)},df);
            sumOfCommonPutALXInputs(i,j) = sum([sum(temp1),sum(temp2)]);
            individualDBXPutALXInputs{sub2ind(size(FiringCorrDBX),i,j)} = [sum(temp1),sum(temp2)];
            
        end
    end
end



