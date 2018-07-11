% Manual Sort List
addpath(genpath('../')); 
ml = readtable('../../SynapseDetector/ManualSortList.csv');
mlMat = table2array(ml);
load('IntConnMatrix.mat')
load('IntPartners.mat')
load('AllCells.mat');
load('ConnMatrixPre.mat');

mlOrdered = [];
mlOrderedTypes = []
mlSize = [];
for i = 1:size(mlMat,2)
    temp = mlMat(:,i)
    temp = temp(~isnan(temp));
    mlSize = [mlSize,size(temp,1)];
    mlOrdered = [mlOrdered; temp];
    mlOrderedTypes = [mlOrderedTypes; repmat(ml.Properties.VariableNames(i),size(temp,1),1)];
end
    
% assemble the ordered Matrix
[m,n,v] = intersect(mlOrdered,AllCells,'stable');

for i =1:size(m,1)
    mlConn(i,:) = ConnMatrixPre(v(i),v);
end

figure;
cspy(mlConn,'Colormap',colorcet('L8'),'Levels',255,'MarkerSize',25);
hold on;
blocks =0;
for i = 1:length(mlSize)
    blocks = blocks+mlSize(i)
    line([0,size(mlConn,1)],[blocks,blocks],'color','k');
    line([blocks,blocks],[0,size(mlConn,1)],'color','k');
end

[a,b] = ismember(functionalCellIDs_new(CellDisplayOrder),AllCells);
b(2) =24; % just a place holder until the cell is validated.
mlInt = ConnMatrixPre(b,v);
figure;
cspy(mlInt,'Colormap',colorcet('L8'),'Levels',255,'MarkerSize',25);
hold on;
%mlSize = [0,mlSize];
blocks = 0;
for i = 1:length(mlSize)
    blocks = blocks+mlSize(i)
   % line([0,size(mlInt,1)],[blocks,blocks],'color','k');
    line([blocks,blocks],[0,size(mlMat,1)],'color','k');
    
end
