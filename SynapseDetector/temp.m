%
clear
df = readtable('/Users/ashwin/Documents/SynapseDetector/04152019.csv');

% load integrator neurons
load('LeadLikeDBX.mat');
load('LagLikeDBX.mat');
load('LeadLikeALX.mat');
load('LagLikeALX.mat');
load('leadNeurons.mat');
load('lagNeurons.mat');


% load Motor neurons
load('smallABDneurons.mat');
load('largeABDneurons.mat');
load('ABDr.mat');
load('ABDc.mat');
load('ABDIr.mat');
load('ABDIc.mat');


% load Saccadic neurons

load('leadSaccadicAxons.mat');
load('lagSaccadicAxons.mat');

leadSaccadicAxons = unique(leadSaccadicAxons);
lagSaccadicAxons  = unique(lagSaccadicAxons);


% load connectome
load('AllCells.mat') ;
load('ConnMatrixPre.mat')

%

confirmedALX = [76181 76201 76187 76184 76192 76197 ];
putativeALX = [76657 76623 76701 76661 76673 76690 76677 76624 77327 77580 ...
    77341 77344 77868 77872 77349 76634 77592 77578 77607 78629 77581 77591 ...
    77602 77821 77844 77845 77822 78406 78667 79022 79033 78404 78441 78421 ...
    79044 79046 79048 80221 78853 79017 79852 78451 79042 80596 80606 78911 ...
    79746 80271 79720 79976 77586 77369 78633 80750 77142 79060 78453 80885];

allALX = [confirmedALX,putativeALX];

confirmedDBX = [76199 76200 76182 76183 76185 76186 76189 76191 76188 ];
putativeDBX = [77582 77588 77595 77597 77598 77599 77605 79040 76399 76828 ...
    76829 78435 76826 76874 76289 76542 76832 76836 76838 76877 76883 76884 ...
    78400 78429 76383 76867 77594 76830 76876 78440 78438 78447 79045 78434 ...
    78450 76878 77683 77450 78448 ];

allDBX = [confirmedDBX, putativeDBX];


confirmedBarhl = [76198, 76190, 76193 ,76194, 76195, 76196];
putativeBarhl = [78452, 80224, 78391];


allBARHL = [confirmedBarhl,putativeBarhl];


 vestibularCellIds = [76631,76700,76656,76655,76660,76679,77393,77395,77375,...
     77815,77803,77807,80591,80579,77255,77697,77364,80223,81126,80760,81117,...
     80672,80622,80572,80569,81122,78426];
 
 MVNs = [76664 76784 76783 77394 77326 77772 77766 77756 77753 77767 77775 ...
     77780 80883 80247 80669 80698 80762 80748 80696 80761 80749 81070 81044 ...
     80991 80967 81008 80883 81045 80986 81023 78647 78619 80247 ];
 
 allVest = [vestibularCellIds,MVNs];
 


 ABDr_CellIDs = [77648, 77710, 77300, 77705, 77305, 77301, 77709, 77672, 77302];
    ABDc_CellIDs = [77154, 77646, 77682 ,77628 ,77295 , 77652 ,77292 ,77688 ,77654 ,77658 ,77657 ,77662, 77296];
    ABDIr_CellIDs = [77631, 77150, 77618, 77886, 78547, 77158, 78556, 78552, 77665, 77668, 77634, 78553];
    ABDIc_CellIDs = [77148, 77625, 77641, 77692, 77144, 77643, 77640, 79051, 79066, 78574];
    
    allABD = [ABDr_CellIDs,ABDc_CellIDs,ABDIr_CellIDs,ABDIc_CellIDs];

%%

cellOrder = [allALX,allDBX,allBARHL,leadSaccadicAxons',lagSaccadicAxons',allVest,allABD];

for i = 1:size(cellOrder,2)
    index(i) = find(cellOrder(i) == AllCells);
end
    
ConnMat = zeros(size(index,2));

for i = 1:size(index,2)
    ConnMat(i,:) = ConnMatrixPre(index(i),index);
end

%%

map = colorcet('L8','N',max(ConnMat(:)));

cspy(ConnMat,'Colormap',map,'Levels',max(ConnMat(:)),'MarkerSize',25);

endALX = find(allALX(end) == AllCells(index));
endDBX = find(allDBX(end) == AllCells(index));
endBARHL = find(allBARHL(end) == AllCells(index));
endLead = find(leadSaccadicAxons(end) == AllCells(index));
endLag = find(lagSaccadicAxons(end) == AllCells(index));
endVest = find(vestibularCellIds(end) == AllCells(index));
endMVN = find(MVNs(end) == AllCells(index));
endABD = find(ABDc_CellIDs(end) == AllCells(index));


hold on;
line([0,size(ConnMat,1)],[endALX,endALX],'color','k');
line([0,size(ConnMat,1)],[endDBX,endDBX],'color','k');
line([0,size(ConnMat,1)],[endBARHL,endBARHL],'color','k');
line([0,size(ConnMat,1)],[endLead,endLead],'color','k');
line([0,size(ConnMat,1)],[endLag,endLag],'color','k');
line([0,size(ConnMat,1)],[endVest,endVest],'color','k');
line([0,size(ConnMat,1)],[endMVN(2),endMVN(2)],'color','k');
line([0,size(ConnMat,1)],[endABD,endABD],'color','k');


line([endALX,endALX],[0,size(ConnMat,1)],'color','k');
line([endDBX,endDBX],[0,size(ConnMat,1)],'color','k');
line([endBARHL,endBARHL],[0,size(ConnMat,1)],'color','k');
line([endLead,endLead],[0,size(ConnMat,1)],'color','k');
line([endLag,endLag],[0,size(ConnMat,1)],'color','k');
line([endVest,endVest],[0,size(ConnMat,1)],'color','k');
line([endMVN(2),endMVN(2)],[0,size(ConnMat,1)],'color','k');
line([endABD,endABD],[0,size(ConnMat,1)],'color','k');

%set(gca,'XTick',[],'YTick',[]);
set(gca,'XTick',1:size(ConnMat,1),'XTickLabel',AllCells(index));

box on;