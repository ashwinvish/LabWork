%DBX analysis
addpath(genpath('../'));
colors = cbrewer('qual', 'Set1', 10);
startup;

% DBX - red
% ALX - blue
% Barhl - green

if ismac
    addpath(genpath('/Users/admin/Documents/Scripts'));
    df = readtable('/Users/ashwin/Documents/SynapseDetector/09202018.csv');
else
    addpath(genpath('/usr/people/ashwinv/seungmount/research/Ashwin/Scripts'));
    df = readtable('/usr/people/ashwinv/seungmount/research/Ashwin/SynapseDetector/09202018.csv');
end

%load('FiringRates.mat');
load('TAU.mat');
load('STA.mat');
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

for i =1:length(DbxCells)
    [DBXpreSynaticPartners{i},DBXpreSynapsePartnerID{i}] = SynapticPartners(DbxCells(i),1,df);
end

for i = 1:length(AlxCells)
    [ALXpreSynaticPartners{i},ALXpreSynapsePartnerID{i}] = SynapticPartners(AlxCells(i),1,df);
end

for i = 1:length(BarhlCells)
    [BARHLpreSynaticPartners{i},BARHLpreSynapsePartnerID{i}] = SynapticPartners(BarhlCells(i),1,df);
end

subplot(4,4,5)
plot(Firing(:,OriginalCellOrderDBX),'color',colors(1,:));
hold on
plot(Firing(:,OriginalCellOrderALX),'color',colors(2,:));
plot(Firing(:,OriginalCellOrderBARHL),'color',colors(3,:));
ylabel('STA df/f');
xlabel('time');
set(gca, 'XColor','w','YColor','w','Color','none');
box off;
axis square;
offsetAxes;



FiringCorrDBX =  corr(Firing(:,OriginalCellOrderDBX));
FiringCorrDBX(1:size(FiringCorrDBX,1)+1:end) = NaN*zeros(1,size(FiringCorrDBX,1));

FiringCorrALX =  corr(Firing(:,OriginalCellOrderALX));
FiringCorrALX(1:size(FiringCorrALX,1)+1:end) = NaN*zeros(1,size(FiringCorrALX,1));

FiringCorrBARHL =  corr(Firing(:,OriginalCellOrderALX));
FiringCorrBARHL(1:size(FiringCorrBARHL,1)+1:end) = NaN*zeros(1,size(FiringCorrBARHL,1));


% Number of common inputs for all pairs

for i = 1:length(DbxCells)
    for j = 1:length(DbxCells)
        CommonInputsDBX{sub2ind(size(FiringCorrDBX),i,j)} = intersect(DBXpreSynaticPartners{i},DBXpreSynaticPartners{j});
    end
end

for i = 1:length(AlxCells)
    for j = 1:length(AlxCells)
        CommonInputsALX{sub2ind(size(FiringCorrALX),i,j)} = intersect(ALXpreSynaticPartners{i},ALXpreSynaticPartners{j});
    end
end

for i = 1:length(BarhlCells)
    for j = 1:length(BarhlCells)
        CommonInputsBARHL{sub2ind(size(FiringCorrBARHL),i,j)} = intersect(BARHLpreSynaticPartners{i},BARHLpreSynaticPartners{j});
    end
end

% Number of Common Inputs

%DBX
for i = 1:size(CommonInputsDBX,2)
    [row,col] = ind2sub(size(FiringCorrDBX),i);
    if row == col
        NumberOfCommonInputsDBX(row,col) = NaN;
    else
        NumberOfCommonInputsDBX(row,col) =  length(CommonInputsDBX{i}) ;
    end
end

%ALX
for i = 1:size(CommonInputsALX,2)
    [row,col] = ind2sub(size(FiringCorrALX),i);
    if row == col
        NumberOfCommonInputsALX(row,col) = NaN;
    else
        NumberOfCommonInputsALX(row,col) =  length(CommonInputsALX{i}) ;
    end
end

%BARHL
for i = 1:size(CommonInputsBARHL,2)
    [row,col] = ind2sub(size(FiringCorrBARHL),i);
    if row == col
        NumberOfCommonInputsBARHL(row,col) = NaN;
    else
        NumberOfCommonInputsBARHL(row,col) =  length(CommonInputsBARHL{i}) ;
    end
end


% plot Number of common inputs vs. Firing coorelations

subplot(4,4,6)
scatter(NumberOfCommonInputsDBX(:),FiringCorrDBX(:),50,'MarkerFaceColor',colors(1,:),'MarkerEdgeColor','none');
hold on;
scatter(NumberOfCommonInputsALX(:),FiringCorrALX(:),50,'MarkerFaceColor',colors(2,:),'MarkerEdgeColor','none');
scatter(NumberOfCommonInputsBARHL(:),FiringCorrBARHL(:),50,'MarkerFaceColor',colors(3,:),'MarkerEdgeColor','none');
ylabel('STA correlation');
xlabel('Number of common inputs');
set(gca, 'XColor','w', 'YColor','w','Color','none');
axis square;
offsetAxes;

%% extract size of presynaptic PSD from common inputs


for i = 1:length(DbxCells)
    for j = 1:length(DbxCells)
        if i==j
            sumOfCommonDBXInputSize(i,j) = NaN;
        else
            [temp1,temp2]= commonInputSize(DbxCells(i),DbxCells(j),CommonInputsDBX{sub2ind(size(FiringCorrDBX),i,j)},df);
            sumOfCommonDBXInputSize(i,j) = sum([sum(temp1),sum(temp2)]);
            clear {'temp1', 'temp2'};
        end
    end
end

for i = 1:length(AlxCells)
    for j = 1:length(AlxCells)
        if i==j
            sumOfCommonALXInputSize(i,j) = NaN;
        else
            [temp1,temp2]= commonInputSize(AlxCells(i),AlxCells(j),CommonInputsALX{sub2ind(size(FiringCorrALX),i,j)},df);
            sumOfCommonALXInputSize(i,j) = sum([sum(temp1),sum(temp2)]);
            clear {'temp1', 'temp2'};
        end
    end
end

for i = 1:length(BarhlCells)
    for j = 1:length(BarhlCells)
        if i==j
            sumOfCommonBARHLInputSize(i,j) = NaN;
        else
            [temp1,temp2]= commonInputSize(BarhlCells(i),BarhlCells(j),CommonInputsBARHL{sub2ind(size(FiringCorrBARHL),i,j)},df);
            sumOfCommonBARHLInputSize(i,j) = sum([sum(temp1),sum(temp2)]);
            clear {'temp1', 'temp2'};        
        end
    end
end

%%

subplot(4,4,7)
scatter(sumOfCommonDBXInputSize(:), FiringCorrDBX(:),50,'MarkerFaceColor',colors(1,:),'MarkerEdgeColor','none');
hold on;
scatter(sumOfCommonALXInputSize(:), FiringCorrALX(:),50,'MarkerFaceColor',colors(2,:),'MarkerEdgeColor','none');
scatter(sumOfCommonBARHLInputSize(:), FiringCorrBARHL(:),50,'MarkerFaceColor',colors(3,:),'MarkerEdgeColor','none');
ylabel('STA correlation');
xlabel('Sum of common Input PSD voxels');
set(gca, 'XColor','w','YColor','w','Color','none');
axis square;
offsetAxes;

%% extract common inputs that are saccadic, vestibular or rest [DBX only]
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

%% extract common inputs that are saccadic, vestibular or rest [ALX only]
for i = 1:size(CommonInputsALX,2)
    [a,b,c] = intersect(mlOrdered,CommonInputsALX{i});
    commonInputsCuratedALX{i} = a;
    commonInputsCuratedTypeALX{i} = mlOrderedTypes(b);
end

% VestibularNeurons
for i =1:size(CommonInputsALX,2)
   index =  find(strcmp(commonInputsCuratedTypeALX{i},'Vestibular'));
   commonALXVestibularInputs{i} = commonInputsCuratedALX{i}(index);  
   clear index;
   index =  find(strcmp(commonInputsCuratedTypeALX{i},'Bushy'));
   commonALXBushyInputs{i} = commonInputsCuratedALX{i}(index);  
   clear index;
   index =  find(strcmp(commonInputsCuratedTypeALX{i},'Lateral_Dorsal'));
   commonALXLateralDInputs{i} = commonInputsCuratedALX{i}(index);
   clear index;
   index =  find(strcmp(commonInputsCuratedTypeALX{i},'PutativeALX'));
   commonALXPutALXInputs{i} = commonInputsCuratedALX{i}(index);
end

for i = 1:length(AlxCells)
    for j = 1:length(AlxCells)
        if i==j
            sumOfCommonVestibularInputs(i,j) = NaN;
            sumOfCommonBushyInputs(i,j) = NaN;
            sumOfCommonLateralDInputs(i,j) = NaN;
            sumOfCommonPutALXInputs(i,j) = NaN;
            individualALXVestibularInputs{sub2ind(size(FiringCorrALX),i,j)} = [NaN, NaN];
            individualALXBushyInputs{sub2ind(size(FiringCorrALX),i,j)} = [NaN, NaN];
            individualALXLateralDInputs{sub2ind(size(FiringCorrALX),i,j)} = [NaN, NaN];
            individualALXPutALXInputs{sub2ind(size(FiringCorrALX),i,j)} = [NaN, NaN];
        else
            [temp1,temp2]= commonInputSize(AlxCells(i),AlxCells(j),commonALXVestibularInputs{sub2ind(size(FiringCorrALX),i,j)},df);
            sumOfCommonVestibularInputs(i,j) = sum([sum(temp1),sum(temp2)]);
            individualALXVestibularInputs{sub2ind(size(FiringCorrALX),i,j)} = [sum(temp1),sum(temp2)];
            clear {temp1,temp2};
            [temp1,temp2]= commonInputSize(AlxCells(i),AlxCells(j),commonALXBushyInputs{sub2ind(size(FiringCorrALX),i,j)},df);
            sumOfCommonBushyInputs(i,j) = sum([sum(temp1),sum(temp2)]);
            individualALXBushyInputs{sub2ind(size(FiringCorrALX),i,j)} = [sum(temp1),sum(temp2)];
            clear {temp1,temp2};
            [temp1,temp2]= commonInputSize(AlxCells(i),AlxCells(j),commonALXLateralDInputs{sub2ind(size(FiringCorrALX),i,j)},df);
            sumOfCommonLateralDInputs(i,j) = sum([sum(temp1),sum(temp2)]);
            individualALXLateralDInputs{sub2ind(size(FiringCorrALX),i,j)} = [sum(temp1),sum(temp2)];
            clear {temp1,temp2};
            [temp1,temp2]= commonInputSize(AlxCells(i),AlxCells(j),commonALXPutALXInputs{sub2ind(size(FiringCorrALX),i,j)},df);
            sumOfCommonPutALXInputs(i,j) = sum([sum(temp1),sum(temp2)]);
            individualALXPutALXInputs{sub2ind(size(FiringCorrALX),i,j)} = [sum(temp1),sum(temp2)];
            
        end
    end
end


%% extract common inputs that are saccadic, vestibular or rest [BARHL only]

for i = 1:size(CommonInputsBARHL,2)
    [a,b,c] = intersect(mlOrdered,CommonInputsBARHL{i});
    commonInputsCuratedBARHL{i} = a;
    commonInputsCuratedTypeBARHL{i} = mlOrderedTypes(b);
end

% VestibularNeurons
for i =1:size(CommonInputsBARHL,2)
   index =  find(strcmp(commonInputsCuratedTypeBARHL{i},'Vestibular'));
   commonBARHLVestibularInputs{i} = commonInputsCuratedBARHL{i}(index);  
   clear index;
   index =  find(strcmp(commonInputsCuratedTypeBARHL{i},'Bushy'));
   commonBARHLBushyInputs{i} = commonInputsCuratedBARHL{i}(index);  
   clear index;
   index =  find(strcmp(commonInputsCuratedTypeBARHL{i},'Lateral_Dorsal'));
   commonBARHLLateralDInputs{i} = commonInputsCuratedBARHL{i}(index);
   clear index;
   index =  find(strcmp(commonInputsCuratedTypeBARHL{i},'PutativeBARHL'));
   commonBARHLPutALXInputs{i} = commonInputsCuratedBARHL{i}(index);
end

for i = 1:length(BarhlCells)
    for j = 1:length(BarhlCells)
        if i==j
            sumOfCommonVestibularInputs(i,j) = NaN;
            sumOfCommonBushyInputs(i,j) = NaN;
            sumOfCommonLateralDInputs(i,j) = NaN;
            sumOfCommonPutBARHLInputs(i,j) = NaN;
            individualBARHLVestibularInputs{sub2ind(size(FiringCorrBARHL),i,j)} = [NaN, NaN];
            individualBARHLBushyInputs{sub2ind(size(FiringCorrBARHL),i,j)} = [NaN, NaN];
            individualBARHLLateralDInputs{sub2ind(size(FiringCorrBARHL),i,j)} = [NaN, NaN];
            individualBARHLPutALXInputs{sub2ind(size(FiringCorrBARHL),i,j)} = [NaN, NaN];
        else
            [temp1,temp2]= commonInputSize(BarhlCells(i),BarhlCells(j),commonBARHLVestibularInputs{sub2ind(size(FiringCorrBARHL),i,j)},df);
            sumOfCommonVestibularInputs(i,j) = sum([sum(temp1),sum(temp2)]);
            individualBARHLVestibularInputs{sub2ind(size(FiringCorrBARHL),i,j)} = [sum(temp1),sum(temp2)];
            clear {temp1,temp2};
            [temp1,temp2]= commonInputSize(BarhlCells(i),BarhlCells(j),commonBARHLBushyInputs{sub2ind(size(FiringCorrBARHL),i,j)},df);
            sumOfCommonBushyInputs(i,j) = sum([sum(temp1),sum(temp2)]);
            individualBARHLBushyInputs{sub2ind(size(FiringCorrBARHL),i,j)} = [sum(temp1),sum(temp2)];
            clear {temp1,temp2};
            [temp1,temp2]= commonInputSize(BarhlCells(i),BarhlCells(j),commonBARHLLateralDInputs{sub2ind(size(FiringCorrBARHL),i,j)},df);
            sumOfCommonLateralDInputs(i,j) = sum([sum(temp1),sum(temp2)]);
            individualBARHLLateralDInputs{sub2ind(size(FiringCorrBARHL),i,j)} = [sum(temp1),sum(temp2)];
            clear {temp1,temp2};
            [temp1,temp2]= commonInputSize(BarhlCells(i),BarhlCells(j),commonBARHLPutALXInputs{sub2ind(size(FiringCorrBARHL),i,j)},df);
            sumOfCommonPutBARHLInputs(i,j) = sum([sum(temp1),sum(temp2)]);
            individualBARHLPutALXInputs{sub2ind(size(FiringCorrBARHL),i,j)} = [sum(temp1),sum(temp2)];
            
        end
    end
end

%%

subplot(4,4,9)

for i =1:size(individualDBXPutALXInputs,2)
    scatter(individualDBXPutALXInputs{i}(1),FiringCorrDBX(ind2sub(size(FiringCorrDBX),i)),50,'MarkerFaceColor',colors(1,:),'MarkerEdgeColor','none');
    hold on
    scatter(individualDBXPutALXInputs{i}(2),FiringCorrDBX(ind2sub(size(FiringCorrDBX),i)),50,'MarkerFaceColor',colors(1,:),'MarkerEdgeColor','none');    
end

for i = 1:size(individualALXPutALXInputs,2)
    scatter(individualALXPutALXInputs{i}(1),FiringCorrALX(ind2sub(size(FiringCorrDBX),i)),50,'MarkerFaceColor',colors(2,:),'MarkerEdgeColor','none');
    hold on;
    scatter(individualALXPutALXInputs{i}(2),FiringCorrALX(ind2sub(size(FiringCorrDBX),i)),50,'MarkerFaceColor',colors(2,:),'MarkerEdgeColor','none');    
end

for i = 1:size(individualBARHLPutALXInputs,2) 
    scatter(individualBARHLPutALXInputs{i}(1),FiringCorrBARHL(ind2sub(size(FiringCorrDBX),i)),50,'MarkerFaceColor',colors(3,:),'MarkerEdgeColor','none');
    hold on;
    scatter(individualBARHLPutALXInputs{i}(2),FiringCorrBARHL(ind2sub(size(FiringCorrDBX),i)),50,'MarkerFaceColor',colors(3,:),'MarkerEdgeColor','none');    
end
ylabel('STA correlation');
xlabel('Sum of putativeALX inputs');
set(gca,'color','none','XColor','k','YColor','k');
axis square;
offsetAxes;


subplot(4,4,10)

for i =1:size(individualDBXVestibularInputs,2)
    scatter(individualDBXVestibularInputs{i}(1),FiringCorrDBX(ind2sub(size(FiringCorrDBX),i)),50,'MarkerFaceColor',colors(1,:),'MarkerEdgeColor','none');
    hold on
    scatter(individualDBXVestibularInputs{i}(2),FiringCorrDBX(ind2sub(size(FiringCorrDBX),i)),50,'MarkerFaceColor',colors(1,:),'MarkerEdgeColor','none');    
end
for i =1:size(individualALXVestibularInputs,2)
    scatter(individualALXVestibularInputs{i}(1),FiringCorrALX(ind2sub(size(FiringCorrALX),i)),50,'MarkerFaceColor',colors(2,:),'MarkerEdgeColor','none');
    hold on
    scatter(individualALXVestibularInputs{i}(2),FiringCorrALX(ind2sub(size(FiringCorrALX),i)),50,'MarkerFaceColor',colors(2,:),'MarkerEdgeColor','none');    
end
for i =1:size(individualBARHLVestibularInputs,2)
    scatter(individualBARHLVestibularInputs{i}(1),FiringCorrBARHL(ind2sub(size(FiringCorrBARHL),i)),50,'MarkerFaceColor',colors(3,:),'MarkerEdgeColor','none');
    hold on
    scatter(individualBARHLVestibularInputs{i}(2),FiringCorrBARHL(ind2sub(size(FiringCorrBARHL),i)),50,'MarkerFaceColor',colors(3,:),'MarkerEdgeColor','none');    
end
ylabel('STA correlation');
xlabel('Sum of Vestibular inputs');
set(gca,'color','none','XColor','k','YColor','k');
axis square;
offsetAxes;

subplot(4,4,11)

for i =1:size(individualDBXLateralDInputs,2)
    scatter(individualDBXLateralDInputs{i}(1),FiringCorrDBX(ind2sub(size(FiringCorrDBX),i)),50,'MarkerFaceColor',colors(1,:),'MarkerEdgeColor','none');
    hold on
    scatter(individualDBXLateralDInputs{i}(2),FiringCorrDBX(ind2sub(size(FiringCorrDBX),i)),50,'MarkerFaceColor',colors(1,:),'MarkerEdgeColor','none');    
end
for i =1:size(individualALXLateralDInputs,2)
    scatter(individualALXLateralDInputs{i}(1),FiringCorrALX(ind2sub(size(FiringCorrALX),i)),50,'MarkerFaceColor',colors(2,:),'MarkerEdgeColor','none');
    hold on
    scatter(individualALXLateralDInputs{i}(2),FiringCorrALX(ind2sub(size(FiringCorrALX),i)),50,'MarkerFaceColor',colors(2,:),'MarkerEdgeColor','none');    
end
for i =1:size(individualBARHLLateralDInputs,2)
    scatter(individualBARHLLateralDInputs{i}(1),FiringCorrBARHL(ind2sub(size(FiringCorrBARHL),i)),50,'MarkerFaceColor',colors(3,:),'MarkerEdgeColor','none');
    hold on
    scatter(individualBARHLLateralDInputs{i}(2),FiringCorrBARHL(ind2sub(size(FiringCorrBARHL),i)),50,'MarkerFaceColor',colors(3,:),'MarkerEdgeColor','none');    
end
ylabel('STA correlation');
xlabel('Sum of LateralD inputs');
set(gca,'color','none','XColor','k','YColor','k');axis square;
offsetAxes;

subplot(4,4,12)

for i =1:size(individualDBXBushyInputs,2)
    scatter(individualDBXBushyInputs{i}(1),FiringCorrDBX(ind2sub(size(FiringCorrDBX),i)),50,'MarkerFaceColor',colors(1,:),'MarkerEdgeColor','none');
    hold on
    scatter(individualDBXBushyInputs{i}(2),FiringCorrDBX(ind2sub(size(FiringCorrDBX),i)),50,'MarkerFaceColor',colors(1,:),'MarkerEdgeColor','none');    
end
for i =1:size(individualALXBushyInputs,2)
    scatter(individualALXBushyInputs{i}(1),FiringCorrALX(ind2sub(size(FiringCorrALX),i)),50,'MarkerFaceColor',colors(2,:),'MarkerEdgeColor','none');
    hold on
    scatter(individualALXBushyInputs{i}(2),FiringCorrALX(ind2sub(size(FiringCorrALX),i)),50,'MarkerFaceColor',colors(2,:),'MarkerEdgeColor','none');    
end
for i =1:size(individualBARHLBushyInputs,2)
    scatter(individualBARHLBushyInputs{i}(1),FiringCorrBARHL(ind2sub(size(FiringCorrBARHL),i)),50,'MarkerFaceColor',colors(3,:),'MarkerEdgeColor','none');
    hold on
    scatter(individualBARHLBushyInputs{i}(2),FiringCorrBARHL(ind2sub(size(FiringCorrBARHL),i)),50,'MarkerFaceColor',colors(3,:),'MarkerEdgeColor','none');    
end
ylabel('STA correlation');
xlabel('Sum of Bushy Saccadic inputs');
set(gca,'color','none','XColor','k','YColor','k','XLim',[0,20000]);axis square;
offsetAxes;
