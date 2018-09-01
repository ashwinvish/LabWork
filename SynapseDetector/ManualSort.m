% Manual Sort List
addpath(genpath('../')); 
ml = readtable('../..//SynapseDetector/ManualSortList-082718.csv');
mlMat = table2array(ml);
load('IntConnMatrix.mat')
load('IntPartners.mat')
load('AllCells.mat');
load('ConnMatrixPre.mat');

% remove cells that are not in AllCells
k = setdiff(mlMat, AllCells)
k = k(~isnan(k));

[a1,a2,a3] = intersect(k,mlMat)
mlMat(a3) = NaN;

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
blocks(1) =0;
for i = 1:length(mlSize)
    blocks(i+1) = blocks(i)+mlSize(i);
    line([0,size(mlConn,1)],[blocks(i+1),blocks(i+1)],'color','k');
    line([blocks(i+1),blocks(i+1)],[0,size(mlConn,1)],'color','k');
    text(repmat(-50,1,1), blocks(i), sprintf('%s\\rightarrow',cell2mat(ml.Properties.VariableNames(i))), 'HorizontalAlignment','center', 'FontWeight','bold');
    set(gca,'YTick',1:size(mlConn),'YTickLabel',AllCells(v),'XTick',1:size(mlConn),...
        'XTickLabel',AllCells(v),'XTickLabelRotation',45,'XAxisLocation','top');
end


% integrators and putative integrator

nonCuratedCells =  setdiff(AllCells, mlOrdered);
nonIntersectingCells = setdiff(1:size(AllCells,1),v);
mlConnOrder = [v;nonIntersectingCells'];
for i = 1:size(ConnMatrixPre,1)    
mlConnAll(i,:) = ConnMatrixPre(mlConnOrder(i),mlConnOrder);
end

figure;
cspy(mlConnAll,'Colormap',colorcet('L8'),'Levels',255,'MarkerSize',25);
hold on;
set(gca,'XTick',1:size(AllCells,1),'XTickLabel',AllCells(mlConnOrder),'XTickLabelRotation',45,'XAxisLocation','top', ...
    'YTick',1:size(AllCells,1),'YTickLabel',AllCells(mlConnOrder),'YTickLabelRotation',45);
line([0,size(ConnMatrixPre,1)],[size(v,1),size(v,1)],'color','k');
line([size(v,1),size(v,1)],[0,size(ConnMatrixPre,1)],'color','k');

box on;

%% arrange connectivity by |Ipsi Contra|
%                          |Ipsi Contra |

contraCells = [];
ipsiCells = [];

for i = 1: size(AllCells,1)
    if isContra(AllCells(i)) ==1 
        contraCells = [contraCells, AllCells(i)];
    else
        ipsiCells = [ipsiCells,AllCells(i)];
    end
end

[~,locIpsi] = intersect(AllCells, ipsiCells); % location of Ipsi cells in AllCells
[~,locCont] = intersect(AllCells,contraCells);% location of Contrs cells in AllCells


[a,b,c] = intersect(mlOrdered,ipsiCells,'stable'); % location of mlOrdered cell in Ipsi cells
[~,d] = setdiff(ipsiCells,ipsiCells(c))

ipsiContraOrder = [c;d;locCont];
%ipsiContraOrder = [locIpsi;locCont];

for i = 1:size(AllCells,1)
    ipsiContraConn(i,:) = ConnMatrixPre(ipsiContraOrder(i),ipsiContraOrder);
end

figure;
cspy(ipsiContraConn,'Colormap',colorcet('L8'),'Levels',255,'MarkerSize',25);
line([0,size(AllCells,1)],[size(locIpsi,1),size(locIpsi,1)],'color','k');
line([size(locIpsi,1),size(locIpsi,1)],[0,size(AllCells,1)],'color','k');
set(gca,'XTick',1:size(AllCells,1),'XTickLabel',AllCells(ipsiContraOrder),'XTickLabelRotation',45,'XAxisLocation','top', ...
    'YTick',1:size(AllCells,1),'YTickLabel',AllCells(ipsiContraOrder),'YTickLabelRotation',45);    


% [a,b] = ismember(functionalCellIDs_new(CellDisplayOrder),AllCells);
% b(2) =24; % just a place holder until the cell is validated.
% mlInt = ConnMatrixPre(b,v);
% figure;
% cspy(mlInt,'Colormap',colorcet('L8'),'Levels',255,'MarkerSize',25);
% hold on;
% %mlSize = [0,mlSize];
% blocks = 0;
% for i = 1:length(mlSize)
%     blocks = blocks+mlSize(i)
%    % line([0,size(mlInt,1)],[blocks,blocks],'color','k');
%     line([blocks,blocks],[0,size(mlMat,1)],'color','k');
%     
% end
