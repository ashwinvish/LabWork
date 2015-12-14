%Lm Axon collateral Pathlenghts

% fish3003_ch1-eported-000.swc

tree1col{1} = [39,32]
tree1col{2} =[18, 28]
tree1col{3} =[18,27]
tree1col{4} =[18,31]
tree1col{5} =[17,30]

for i = 1:numel(tree1col)
    axonNodes = AxonQueryNodes(tree{1},tree1col{i});
    colPathLengthTemp{i} = findPathLength('fish3003_ch1-exported-000.swc',tree{1},[0.36,0.36,2],axonNodes);
    tree1colPathLength{i} = sum(diff(sort(colPathLengthTemp{i})));
end

clear axonNodes;
clear colPathLengthTemp;
%%
% fish2017_ch2.swc

tree2col{1} = [11,22]
tree2col{2} = [12,21]
tree2col{3} = [12,20]
tree2col{4} = [10,23]

for i = 1:numel(tree2col)
    axonNodes = AxonQueryNodes(tree{2},tree2col{i});
    colPathLengthTemp{i} = findPathLength('fish2017_ch2.swc',tree{2},[0.36,0.36,2],axonNodes);
    tree2colPathLength{i} = sum(diff(sort(colPathLengthTemp{i})));
end

clear axonNodes;
clear colPathLengthTemp;
%%
% fish1059_ch2-axons.swc
 
tree3col{1} = [4,5]
tree3col{2} = [4,12]
tree3col{3} = [3,13]
tree3col{4} = [6,18]
tree3col{5} = [7, 17]
tree3col{6} = [8,16]
tree3col{7} = [9,15]
tree3col{8} = [10,14]
tree3col{9} = [10,11]

for i = 1:numel(tree3col)
    axonNodes = AxonQueryNodes(tree{3},tree3col{i});
    colPathLengthTemp{i} = findPathLength('fish1059_ch2-axons.swc',tree{3},[0.36,0.36,2],axonNodes);
    tree3colPathLength{i} = sum(diff(sort(colPathLengthTemp{i})));
end

clear axonNodes;
clear colPathLengthTemp;

%%
tree4colPathLength = 0;

%%
%fish1013_ch2.swc

tree5col{1} = [31, 54]
tree5col{2} =[30,53]
tree5col{3} =[29,52]
tree5col{4} =[28,47]
tree5col{5} =[27,51]
tree5col{6} =[26,50]
tree5col{7} =[25,19]
tree5col{8} =[48,63]
tree5col{9} =[11,12]
tree5col{10} =[33,55]
tree5col{11} =[34,58]
tree5col{12} =[35,36]
tree5col{13} =[56,57]

for i = 1:numel(tree5col)
    axonNodes = AxonQueryNodes(tree{5},tree5col{i});
    colPathLengthTemp{i} = findPathLength('fish1013_ch2.swc',tree{5},[0.36,0.36,2],axonNodes);
    tree5colPathLength{i} = sum(diff(sort(colPathLengthTemp{i})));
end

clear axonNodes;
clear colPathLengthTemp;

%%

meanLMAxonColPathLengths = [mean(cell2mat(tree1colPathLength)), mean(cell2mat(tree2colPathLength)), mean(cell2mat(tree3colPathLength)), 0, mean(cell2mat(tree5colPathLength))];
SDLMAxonColPathLengths = [std(cell2mat(tree1colPathLength)), std(cell2mat(tree2colPathLength)), std(cell2mat(tree3colPathLength)),0, std(cell2mat(tree5colPathLength))];

x = 1:5;
hold on;

plot(repmat(x(1),1,size(cell2mat(tree1colPathLength),2)), cell2mat(tree1colPathLength), 'Marker','o','MarkerSize',20 ,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerEdgeColor','k', 'LineStyle','none' );
plot(repmat(x(2),1,size(cell2mat(tree2colPathLength),2)), cell2mat(tree2colPathLength), 'Marker','o','MarkerSize',20 ,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerEdgeColor','k', 'LineStyle','none' );
plot(repmat(x(3),1,size(cell2mat(tree3colPathLength),2)), cell2mat(tree3colPathLength), 'Marker','o','MarkerSize',20 ,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerEdgeColor','k', 'LineStyle','none' );
plot(x(4),tree4colPathLength, 'Marker','o','MarkerSize',20 ,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerEdgeColor','k', 'LineStyle','none' );
plot(repmat(x(5),1,size(cell2mat(tree5colPathLength),2)), cell2mat(tree5colPathLength), 'Marker','o','MarkerSize',20 ,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerEdgeColor','k', 'LineStyle','none' );

plot(1:5 , meanLMAxonColPathLengths, 'Marker','o','MarkerSize',20 ,'MarkerFaceColor','r','MarkerEdgeColor','k', 'LineStyle','none' );
plot([1:5;1:5], [(meanLMAxonColPathLengths-SDLMAxonColPathLengths) ; (meanLMAxonColPathLengths+SDLMAxonColPathLengths)], 'Color','k','LineWidth',2);

set(gca,'XLim', [1 5] , 'XTick', 1:6, 'FontName', 'Arial', 'FontSize', 20);
set(gcf,'color','w');
xlabel('Neuron #', 'FontName', 'Arial', 'FontSize', 20);
ylabel('Collateral Path Length in \mum', 'FontName', 'Arial', 'FontSize', 20);
box off;
axis square;
hold off;

