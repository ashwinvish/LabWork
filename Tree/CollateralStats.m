% Axon Collateral Stats

% Lm Collaterals
figure();
meanLMAxonColPathLengths = [mean(cell2mat(tree1colPathLength)), mean(cell2mat(tree2colPathLength)), mean(cell2mat(tree3colPathLength)), 0, mean(cell2mat(tree5colPathLength))];
SDLMAxonColPathLengths = [std(cell2mat(tree1colPathLength)), std(cell2mat(tree2colPathLength)), std(cell2mat(tree3colPathLength)),0, std(cell2mat(tree5colPathLength))];

x = 1:5;
hold on;

plot(repmat(x(1),1,size(cell2mat(tree1colPathLength),2)), cell2mat(tree1colPathLength), 'Marker','o','MarkerSize',25 ,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerEdgeColor','k', 'LineStyle','none' );
plot(repmat(x(2),1,size(cell2mat(tree2colPathLength),2)), cell2mat(tree2colPathLength), 'Marker','o','MarkerSize',25 ,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerEdgeColor','k', 'LineStyle','none' );
plot(repmat(x(3),1,size(cell2mat(tree3colPathLength),2)), cell2mat(tree3colPathLength), 'Marker','o','MarkerSize',25 ,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerEdgeColor','k', 'LineStyle','none' );
plot(x(4),tree4colPathLength, 'Marker','o','MarkerSize',25 ,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerEdgeColor','k', 'LineStyle','none' );
plot(repmat(x(5),1,size(cell2mat(tree5colPathLength),2)), cell2mat(tree5colPathLength), 'Marker','o','MarkerSize',25 ,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerEdgeColor','k', 'LineStyle','none' );

plot(1:5 , meanLMAxonColPathLengths, 'Marker','o','MarkerSize',35 ,'MarkerFaceColor','r','MarkerEdgeColor','k', 'LineStyle','none' );
plot([1:5;1:5], [(meanLMAxonColPathLengths-SDLMAxonColPathLengths) ; (meanLMAxonColPathLengths+SDLMAxonColPathLengths)], 'Color','k','LineWidth',2);

set(gca,'XLim', [1 5] , 'XTick', 1:6, 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
set(gcf,'color','w');
xlabel('Neuron #', 'FontName', 'Arial', 'FontSize', 40);
ylabel('LM collateral pathlength in \mum', 'FontName', 'Arial', 'FontSize', 40);
box off;
axis square;
hold off;

% EM Axon collateral path lengths
figure();
x = 1:7;
hold on;

plot(repmat(x(1),1,size(cell2mat(Int1_4colPathLength),2)),cell2mat(Int1_4colPathLength)/1000,'Marker','o','MarkerSize',25 ,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerEdgeColor','k', 'LineStyle','none' );
plot(repmat(x(2),1,size(cell2mat(Int1_5colPathLength),2)),cell2mat(Int1_5colPathLength)/1000,'Marker','o','MarkerSize',25 ,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerEdgeColor','k', 'LineStyle','none' );
plot(repmat(x(3),1,size(cell2mat(Int1_6colPathLength),2)),cell2mat(Int1_6colPathLength)/1000,'Marker','o','MarkerSize',25 ,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerEdgeColor','k' , 'LineStyle','none' );
plot(repmat(x(4),1,size(cell2mat(Int1_7colPathLength),2)),cell2mat(Int1_7colPathLength)/1000,'Marker','o','MarkerSize',25 ,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerEdgeColor','k' , 'LineStyle','none' );
plot(repmat(x(5),1,size(cell2mat(Int2_9colPathLength),2)),cell2mat(Int2_9colPathLength)/1000,'Marker','o','MarkerSize',25 ,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerEdgeColor','k' , 'LineStyle','none' );
plot(repmat(x(6),1,size(cell2mat(Int3_5colPathLength),2)),cell2mat(Int3_5colPathLength)/1000,'Marker','o','MarkerSize',25 ,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerEdgeColor','k' , 'LineStyle','none' );

meanColPathLengths = [mean(cell2mat(Int1_4colPathLength)),mean(cell2mat(Int1_5colPathLength)), mean(cell2mat(Int1_6colPathLength)),mean(cell2mat(Int1_7colPathLength)),mean(cell2mat(Int2_9colPathLength)),mean(cell2mat(Int3_5colPathLength))];
SDColPathLengths = [std(cell2mat(Int1_4colPathLength)),std(cell2mat(Int1_5colPathLength)), std(cell2mat(Int1_6colPathLength)),std(cell2mat(Int1_7colPathLength)),std(cell2mat(Int2_9colPathLength)),std(cell2mat(Int3_5colPathLength))];
plot(1:6 , meanColPathLengths/1000, 'Marker','o','MarkerSize',35 ,'MarkerFaceColor','r','MarkerEdgeColor','k', 'LineStyle','none' );
plot([1:6;1:6], [(meanColPathLengths-SDColPathLengths)/1000 ; (meanColPathLengths+SDColPathLengths)/1000], 'Color','k','LineWidth',2);  

set(gca,'XLim', [1 6] , 'XTick', 1:6, 'YLim',[0 100], 'FontName', 'Arial', 'FontSize', 40,'LineWidth',2);
set(gcf,'color','w');
xlabel('Neuron #', 'FontName', 'Arial', 'FontSize', 40);
ylabel('EM collateral pathlength in \mum', 'FontName', 'Arial', 'FontSize', 40);
box off;
axis square;
hold off;

% total collateral lenghts

figure();
hold on;

totalEMAxonColPathLengths = [sum(cell2mat(Int1_4colPathLength)), sum(cell2mat(Int1_5colPathLength)), sum(cell2mat(Int1_6colPathLength)), sum(cell2mat(Int1_7colPathLength)), sum(cell2mat(Int2_9colPathLength)), sum(cell2mat(Int3_5colPathLength))];
totalLMAxonColPathLengths = [sum(cell2mat(tree1colPathLength)), sum(cell2mat(tree2colPathLength)), sum(cell2mat(tree3colPathLength)), 0, sum(cell2mat(tree5colPathLength))];

plot(repmat(0.5,1,size(totalEMAxonColPathLengths,2)),totalEMAxonColPathLengths/1000,'Marker','o','MarkerSize',25 ,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerEdgeColor','k', 'LineStyle','none' );
h = plot(repmat(1,1,size(totalLMAxonColPathLengths,2)),totalLMAxonColPathLengths,'Marker','o','MarkerSize',25 ,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerEdgeColor','k', 'LineStyle','none' )
% drawnow;
% hMarkers = h.MarkerHandle;
% hMarkers.FaceColorData =  uint8(255*[0,1,0,0.1])';  % Alpha=0.5 => 50% transparent red
plot(0.5, mean(totalEMAxonColPathLengths/1000),'Marker','o','MarkerSize',35 ,'MarkerFaceColor','k','MarkerEdgeColor','k', 'LineStyle','none' );
plot(1, mean(totalLMAxonColPathLengths),'Marker','o','MarkerSize',35 ,'MarkerFaceColor','k','MarkerEdgeColor','k', 'LineStyle','none' );
plot([0.5;0.5], [ mean(totalEMAxonColPathLengths/1000) + std(totalEMAxonColPathLengths/1000); mean(totalEMAxonColPathLengths/1000) - std(totalEMAxonColPathLengths/1000)],'Color','k','LineWidth',2);
plot([1;1], [ mean(totalLMAxonColPathLengths) + std(totalLMAxonColPathLengths);mean(totalLMAxonColPathLengths) - std(totalLMAxonColPathLengths)],'Color','k','LineWidth',2);
pbaspect([1 2 1]);
set(gca,'XLim', [0 1.5], 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
set(gcf,'color','w');
set(gca,'XTick', [0.5,1],'XTickLabel', {'EM'; 'LM'}, 'FontName', 'Arial', 'FontSize', 40);
ylabel('Total collateral pathlength in \mum', 'FontName', 'Arial', 'FontSize', 40);
box off;
%axis square;
hold off;


% total axonal Lengths

figure();
hold on;

TotalEMAxonPathLengths = [axLength(4),axLength(5),axLength(6),axLength(7), axLength(16), axLength(21)];
TotalEMAxonPathLengths = TotalEMAxonPathLengths/1000; % covert to um
TotaLMAxonPathLengths = cell2mat(LMAxPathLen);

plot(repmat(0.5,1,size(TotalEMAxonPathLengths,2)), TotalEMAxonPathLengths,'Marker','o','MarkerSize',25 ,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerEdgeColor','k', 'LineStyle','none' );
h1 = plot(repmat(1,1,size(TotaLMAxonPathLengths,2) ),TotaLMAxonPathLengths ,'Marker','o','MarkerSize',25 ,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerEdgeColor','k', 'LineStyle','none');
% drawnow;
% h1Markers = h1.MarkerHandle;
% h1Markers.FaceColorData =  uint8(255*[0,1,0,0.1])';  % Alpha=0.5 => 50% transparent red
plot(0.5, mean(TotalEMAxonPathLengths),'Marker','o','MarkerSize',35 ,'MarkerFaceColor','k','MarkerEdgeColor','k', 'LineStyle','none' );
plot(1, mean(TotaLMAxonPathLengths),'Marker','o','MarkerSize',35 ,'MarkerFaceColor','k','MarkerEdgeColor','k', 'LineStyle','none' );
plot([0.5;0.5], [ mean(TotalEMAxonPathLengths) + std(TotalEMAxonPathLengths); mean(TotalEMAxonPathLengths) - std(TotalEMAxonPathLengths)],'Color','k','LineWidth',2);
plot([1;1], [ mean(TotaLMAxonPathLengths) + std(TotaLMAxonPathLengths);mean(TotaLMAxonPathLengths) - std(TotaLMAxonPathLengths)],'Color','k','LineWidth',2);
pbaspect([1 2 1]);

set(gca,'XLim', [0 1.5], 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
set(gcf,'color','w');
set(gca,'XTick', [0.5,1],'XTickLabel', {'EM'; 'LM'}, 'FontName', 'Arial', 'FontSize', 40);
ylabel('Total axonal pathlength in \mum', 'FontName', 'Arial', 'FontSize', 40);
box off;
%axis square;
hold off;


%number of Collaterals
figure();
hold on;

NumberOfLMCollaterals = [numel(tree1col),numel(tree2col),numel(tree3col),0,numel(tree5col)];
NumberOfEMCollaterals =  [numel(Int1_4col), numel(Int1_5col), numel(Int1_6col), numel(Int1_7col), numel(Int2_9col), numel(Int3_5col)];

plot(0.5*ones(1,size(NumberOfEMCollaterals,2)), NumberOfEMCollaterals, 'Marker','o','MarkerSize',25 ,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerEdgeColor','k', 'LineStyle','none' );
plot(ones(1,size(NumberOfLMCollaterals,2)),NumberOfLMCollaterals , 'Marker','o','MarkerSize',25 ,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerEdgeColor','k', 'LineStyle','none' );
plot(0.5, mean(NumberOfEMCollaterals), 'Marker','o','MarkerSize',35 ,'MarkerFaceColor','k','MarkerEdgeColor','k', 'LineStyle','none' );
plot(1, mean(NumberOfLMCollaterals) , 'Marker','o','MarkerSize',35 ,'MarkerFaceColor','k','MarkerEdgeColor','k', 'LineStyle','none' );
plot([0.5;0.5], [mean(NumberOfEMCollaterals)+std(NumberOfEMCollaterals);mean(NumberOfEMCollaterals)-std(NumberOfEMCollaterals)], 'Color','k','LineWidth',2);
plot([1;1], [mean(NumberOfLMCollaterals)+std(NumberOfLMCollaterals);mean(NumberOfLMCollaterals)-std(NumberOfLMCollaterals)], 'Color','k','LineWidth',2);
pbaspect([1 2 1]);

set(gca,'XLim', [0 1.5], 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
set(gcf,'color','w');
set(gca,'XTick', [0.5,1],'XTickLabel', {'EM'; 'LM'}, 'FontName', 'Arial', 'FontSize', 40);
ylabel('Total number of collaterals', 'FontName', 'Arial', 'FontSize', 40);
box off;
%axis square;
hold off;







