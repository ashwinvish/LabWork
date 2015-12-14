% AxonCollateral lengths for trees

% Int1_5
Int1_5col{1} = [155,154]
Int1_5col{2} = [152,155]
Int1_5col{3} = [152,156,157,158]
Int1_5col{4} = [140,143]
Int1_5col{5} = [132,134]
Int1_5col{6} = [128,131,136,141,142,144]
Int1_5col{7} = [103,117,123,125,126,127,137,139,150,151,153]
Int1_5col{8} = [124,117,131,132,140,152]
Int1_5col{9} = [92,99,100,103,115]
Int1_5col{10} = [77,92,106,110,114,116,133,135,145,146,147,148,149]
Int1_5col{11} = [77,80,81,82,84,85,86]
Int1_5col{12} = [163,165,168]
Int1_5col{13} = [163,165,166,167,169,170]

LR_Caudal = [56,107,159];
LR_RostrAL = [56,77];


for i = 1:numel(Int1_5col)
    axonNodes = AxonQueryNodes(allTrees{5},Int1_5col{i});
    colPathLengthTemp{i} = findPathLength('Int1_5_WithTags.swc',allTrees{5},[5,5,45],axonNodes);
    Int1_5colPathLength{i} = sum(diff(sort(colPathLengthTemp{i})));
end

clear axonNodes;
clear colPathLengthTemp;


%%
% Int1_4
Int1_4col{1} = [67,68,70]
Int1_4col{2} = [62,65,66,68,71,72,73,74]
Int1_4col{3} = [62,63,65,69]
Int1_4col{4} = [48,49,56]
Int1_4col{5} = [48,51]

Int1_4LR_Rostral = [41,48]
Int1_4LR_Caudal = [41,51,52,53,54,55,57,58,59,60,61]

for i = 1:numel(Int1_4col)
    axonNodes = AxonQueryNodes(allTrees{4},Int1_4col{i});
    colPathLengthTemp{i} = findPathLength('Int1_4_WithTags.swc',allTrees{4},[5,5,45],axonNodes);
    Int1_4colPathLength{i} = sum(diff(sort(colPathLengthTemp{i})));
end

clear axonNodes;
clear colPathLengthTemp;
 
%%
% Int1_6

Int1_6col{1} = [130,138]
Int1_6col{2} = [130,137,145]
Int1_6col{3} = [137,147]
Int1_6col{4} = [74,118,130]
Int1_6col{5} = [136,146];
Int1_6col{6} = [74,118,120,121,136]
Int1_6col{7} = [120,136,162]
Int1_6col{8} = [162,168]
Int1_6col{9} = [162,169,170,171,172,173]
Int1_6col{10} = [148,155]
Int1_6col{11} = [148,150,151,157]
Int1_6col{12} = [144,149,152,153]
Int1_6col{13} = [144,134,148]
Int1_6col{14} = [134,139]
Int1_6col{15} = [123,126,128,133]
Int1_6col{16} = [159,165]
Int1_6col{17} = [159,167]
Int1_6col{18} = [160,161]
Int1_6col{19} = [160,164]
Int1_6col{20} = [158,163]
Int1_6col{21} = [131,135,166]
Int1_6col{22} = [135,141,143,154]

Int1_6LR_Caudal = [82,96,98,99,134]
Int1_6LR_Rostral = [131,140,142,156,82,123,129]
Int1_6LR_Caudal2 = [74,118,120,121,136,162]
Int1_6LR_Rostral2 = [74,118,130]

for i = 1:numel(Int1_6col)
    axonNodes = AxonQueryNodes(allTrees{6},Int1_6col{i});
    colPathLengthTemp{i} = findPathLength('Int1_6_WithTags.swc',allTrees{6},[5,5,45],axonNodes);
    Int1_6colPathLength{i} = sum(diff(sort(colPathLengthTemp{i})));
end

clear axonNodes;
clear colPathLengthTemp;

%%
% Int1_7

Int1_7col{1} = [115,129]
Int1_7col{2} = [115,120]
Int1_7col{3} = [81,115,99]
Int1_7col{4} = [30,75,81]
Int1_7col{5} = [75,110,114]
Int1_7col{6} = [110,116,117,121,122,123]

Int1_7LR_Contra = [30,76]

for i = 1:numel(Int1_7col)
    axonNodes = AxonQueryNodes(allTrees{7},Int1_7col{i});
    colPathLengthTemp{i} = findPathLength('Int1_7_WithTags.swc',allTrees{7},[5,5,45],axonNodes);
    Int1_7colPathLength{i} = sum(diff(sort(colPathLengthTemp{i})));
end

clear axonNodes;
clear colPathLengthTemp;

%%
% Int 2_9

Int2_9col{1} = [78,81]
Int2_9col{2} = [43,55]

Int2_9LR_Rostal = [17,47]
Int2_9LR_Caudal = [78,79,80,82,83,84,85,86,87,88,17,66,43,67]


for i = 1:numel(Int2_9col)
    axonNodes = AxonQueryNodes(allTrees{16},Int2_9col{i});
    colPathLengthTemp{i} = findPathLength('Int2_9_WithTags.swc',allTrees{16},[5,5,45],axonNodes);
    Int2_9colPathLength{i} = sum(diff(sort(colPathLengthTemp{i})));
end

clear axonNodes;
clear colPathLengthTemp;

%%
% Int3_5

Int3_5col{1} =  [105,108]
Int3_5col{2} =  [105,106,109]
Int3_5col{3} = [105,106,110]
Int3_5col{4} = [99,100]
Int3_5col{5} = [99,105]
Int3_5col{6} = [87,89]
Int3_5col{7} = [81,84]
Int3_5col{8} = [80,82]
Int3_5col{9} = [83,92]
Int3_5col{10} = [101,102]
Int3_5col{11} = [88,93,90,103,104,107]
 
Int3_5LR_Caudal = [99,96,95,87,81,52,34]
Int3_5LR_Rostral = [52,80,83,86,93,101,102,111,112,113,114,115]



for i = 1:numel(Int3_5col)
    axonNodes = AxonQueryNodes(allTrees{21},Int3_5col{i});
    colPathLengthTemp{i} = findPathLength('Int3_5_WithTags.swc',allTrees{21},[5,5,45],axonNodes);
    Int3_5colPathLength{i} = sum(diff(sort(colPathLengthTemp{i})));
end

clear axonNodes;
clear colPathLengthTemp;


%% plots

x = 1:7;
hold on;

plot(repmat(x(1),1,size(cell2mat(Int1_4colPathLength),2)),cell2mat(Int1_4colPathLength)/1000,'Marker','o','MarkerSize',20 ,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerEdgeColor','k', 'LineStyle','none' );
plot(repmat(x(2),1,size(cell2mat(Int1_5colPathLength),2)),cell2mat(Int1_5colPathLength)/1000,'Marker','o','MarkerSize',20 ,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerEdgeColor','k', 'LineStyle','none' );
plot(repmat(x(3),1,size(cell2mat(Int1_6colPathLength),2)),cell2mat(Int1_6colPathLength)/1000,'Marker','o','MarkerSize',20 ,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerEdgeColor','k' , 'LineStyle','none' );
plot(repmat(x(4),1,size(cell2mat(Int1_7colPathLength),2)),cell2mat(Int1_7colPathLength)/1000,'Marker','o','MarkerSize',20 ,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerEdgeColor','k' , 'LineStyle','none' );
plot(repmat(x(5),1,size(cell2mat(Int2_9colPathLength),2)),cell2mat(Int2_9colPathLength)/1000,'Marker','o','MarkerSize',20 ,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerEdgeColor','k' , 'LineStyle','none' );
plot(repmat(x(6),1,size(cell2mat(Int3_5colPathLength),2)),cell2mat(Int3_5colPathLength)/1000,'Marker','o','MarkerSize',20 ,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerEdgeColor','k' , 'LineStyle','none' );

meanColPathLengths = [mean(cell2mat(Int1_4colPathLength)),mean(cell2mat(Int1_5colPathLength)), mean(cell2mat(Int1_6colPathLength)),mean(cell2mat(Int1_7colPathLength)),mean(cell2mat(Int2_9colPathLength)),mean(cell2mat(Int3_5colPathLength))];
SDColPathLengths = [std(cell2mat(Int1_4colPathLength)),std(cell2mat(Int1_5colPathLength)), std(cell2mat(Int1_6colPathLength)),std(cell2mat(Int1_7colPathLength)),std(cell2mat(Int2_9colPathLength)),std(cell2mat(Int3_5colPathLength))];
plot(1:6 , meanColPathLengths/1000, 'Marker','o','MarkerSize',20 ,'MarkerFaceColor','r','MarkerEdgeColor','k', 'LineStyle','none' );
plot([1:6;1:6], [(meanColPathLengths-SDColPathLengths)/1000 ; (meanColPathLengths+SDColPathLengths)/1000], 'Color','k','LineWidth',2);  

set(gca,'XLim', [1 6] , 'XTick', 1:6, 'YLim',[0 100], 'FontName', 'Arial', 'FontSize', 20);
set(gcf,'color','w');
xlabel('Neuron #', 'FontName', 'Arial', 'FontSize', 20);
ylabel('Collateral Path Length in \mum', 'FontName', 'Arial', 'FontSize', 20);
box off;
axis square;
hold off;
% 