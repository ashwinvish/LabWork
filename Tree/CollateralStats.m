% Axon Collateral Stats

% Lm Collaterals
figure();
meanLMAxonColPathLengths = [mean(cell2mat(tree1colPathLength)), mean(cell2mat(tree2colPathLength)), mean(cell2mat(tree3colPathLength)), mean(cell2mat(tree6colPathLength))];
SDLMAxonColPathLengths = [std(cell2mat(tree1colPathLength)), std(cell2mat(tree2colPathLength)), std(cell2mat(tree3colPathLength)), std(cell2mat(tree6colPathLength))];

x = 1:4;
hold on;

plot(repmat(x(1),1,size(cell2mat(tree1colPathLength),2)), cell2mat(tree1colPathLength), 'Marker','o','MarkerSize',20 ,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerEdgeColor','k', 'LineStyle','none' );
plot(repmat(x(2),1,size(cell2mat(tree2colPathLength),2)), cell2mat(tree2colPathLength), 'Marker','o','MarkerSize',20 ,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerEdgeColor','k', 'LineStyle','none' );
plot(repmat(x(3),1,size(cell2mat(tree3colPathLength),2)), cell2mat(tree3colPathLength), 'Marker','o','MarkerSize',20 ,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerEdgeColor','k', 'LineStyle','none' );
plot(repmat(x(4),1,size(cell2mat(tree6colPathLength),2)), cell2mat(tree6colPathLength), 'Marker','o','MarkerSize',20 ,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerEdgeColor','k', 'LineStyle','none' );

plot(1:4 , meanLMAxonColPathLengths, 'Marker','o','MarkerSize',20 ,'MarkerFaceColor','r','MarkerEdgeColor','k', 'LineStyle','none' );
plot([1:4;1:4], [(meanLMAxonColPathLengths-SDLMAxonColPathLengths) ; (meanLMAxonColPathLengths+SDLMAxonColPathLengths)], 'Color','k','LineWidth',2);

set(gca,'XLim', [1 4] , 'XTick', 1:6, 'FontName', 'Arial', 'FontSize', 20);
set(gcf,'color','w');
xlabel('Neuron #', 'FontName', 'Arial', 'FontSize', 20);
ylabel('LM Collateral Path Length in \mum', 'FontName', 'Arial', 'FontSize', 20);
box off;
axis square;
hold off;

% EM Axon collateral path lengths
figure();
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
ylabel('EM Collateral Path Length in \mum', 'FontName', 'Arial', 'FontSize', 20);
box off;
axis square;
hold off;

% total collateral lenghts

figure();
hold on;

totalEMAxonColPathLengths = [sum(cell2mat(Int1_4colPathLength)), sum(cell2mat(Int1_5colPathLength)), sum(cell2mat(Int1_6colPathLength)), sum(cell2mat(Int1_7colPathLength)), sum(cell2mat(Int2_9colPathLength)), sum(cell2mat(Int3_5colPathLength))];
totalLMAxonColPathLengths = [sum(cell2mat(tree1colPathLength)), sum(cell2mat(tree2colPathLength)), sum(cell2mat(tree3colPathLength)), sum(cell2mat(tree6colPathLength))];

plot(repmat(0.5,1,size(totalEMAxonColPathLengths,2)),totalEMAxonColPathLengths/1000,'Marker','o','MarkerSize',20 ,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerEdgeColor','k', 'LineStyle','none' );
plot(repmat(1,1,size(totalLMAxonColPathLengths,2)),totalLMAxonColPathLengths,'Marker','o','MarkerSize',20 ,'MarkerFaceColor',[0.1,0.75,0.1],'MarkerEdgeColor','k', 'LineStyle','none' );
plot(0.5, mean(totalEMAxonColPathLengths/1000),'Marker','o','MarkerSize',20 ,'MarkerFaceColor','r','MarkerEdgeColor','k', 'LineStyle','none' );
plot(1, mean(totalLMAxonColPathLengths),'Marker','o','MarkerSize',20 ,'MarkerFaceColor','g','MarkerEdgeColor','k', 'LineStyle','none' )
plot([0.5;0.5], [ mean(totalEMAxonColPathLengths/1000) + std(totalEMAxonColPathLengths/1000); mean(totalEMAxonColPathLengths/1000) - std(totalEMAxonColPathLengths/1000)],'Color','k','LineWidth',2);
plot([1;1], [ mean(totalLMAxonColPathLengths) + std(totalLMAxonColPathLengths);mean(totalLMAxonColPathLengths) - std(totalLMAxonColPathLengths)],'Color','k','LineWidth',2);

set(gca,'XLim', [0 1.5], 'FontName', 'Arial', 'FontSize', 20);
set(gcf,'color','w');
set(gca,'XTick', [0.5,1],'XTickLabel', {'EM'; 'LM'}, 'FontName', 'Arial', 'FontSize', 20);
ylabel('Collateral Path Length in \mum', 'FontName', 'Arial', 'FontSize', 20);
box off;
axis square;
hold off;



