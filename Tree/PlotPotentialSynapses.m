%% 

PotSites = [];
PotSitesTemp = [];
index = 1;
cellIDAlxSparse = {'Int1_4','Int1_5','Int1_6','Int1_7','Int2_9','Int3_5','Int3_6'};

for i = 1:numel(cellIDs)
    if ismember(cellIDs(i),cellIDAlxSparse) == 1
        i
        %myCell{index} = [];
        figure();
        
        for j = 1:numel(cellIDs)
            denTree = [];
            subplot(3,8,j);
            axonTree =  UniqueSites(allTrees{i},cellIDs,i,true,false,0);
             pause (1);
            DenTree = allTrees{j};                                        % iterating through all trees
            validNodes =  eval([cellIDs{j},'_axon']);                     % keeping track of axonal nodes
            if isempty(validNodes)                                        % if no axon check
                validNodes = 1:numel(DenTree);
            end                                                           % if no axon then iterate through all nodes
            for jj = 1:numel(DenTree)
                children = DenTree{jj}{2};                                % Consider childern of jj node
                for nn = 1:numel(children)                                % iterate over all children
                    DnTempx = [DenTree{children(nn)}{3}(1); DenTree{children(nn)}{4}{1}(:,1)];
                    DnTempy = [DenTree{children(nn)}{3}(2); DenTree{children(nn)}{4}{1}(:,2)];
                    DnTempz = -[DenTree{children(nn)}{3}(3); DenTree{children(nn)}{4}{1}(:,3)];
                    if length(validNodes)< length(DenTree) && ismember(jj,validNodes)
                        continue;
                    else
                         hold on;
                        h2 = plot3(DnTempx,DnTempy,DnTempz,'-','color',[0.8 0.8 0.8]);
                       
                    end
                    denTree = [denTree; DnTempx DnTempy DnTempz];
                end
            end
            
            [PotSitesTemp] = PotentialSites(axonTree,denTree,i,j);
            
            if ~isempty(PotSitesTemp)
                scatter3(PotSitesTemp(:,1),PotSitesTemp(:,2), PotSitesTemp(:,3), 'Marker','o', 'MarkerFaceColor', 'b');
                str = sprintf('Potential synapses %d',size(PotSitesTemp,1));
                title(str);
            end
            
            myCell{index}{j} = PotSitesTemp;  
            clear PotSitesTemp
            clear denTree;
            clear DenTree;
            clear validNodes;
            clear children;
            
        end
        index = index+1;
    end
end

%%

    for i = 1:numel(myCell)
        myCellSize = cellfun(@size,myCell{1,i},'uni',false);
        for j = 1:numel(myCellSize)
           sizeOfMyCell(i,j) =  myCellSize{1,j}(1);
        end
        nonZeroElements(i) = numel(find(sizeOfMyCell(i,:)));
        sumOfPotentialSynapses(i) = sum(sizeOfMyCell(i,:));
    end
    
   [i,j,s] = find(sizeOfMyCell);
    plot(i,s,'o', 'MarkerSize' , 10, 'MarkerFaceColor','k', 'MarkerEdgeColor', 'none');
    hold on;
    meanNoJitter = sumOfPotentialSynapses./nonZeroElements;
    meanNoJitter(1) = 0;
    plot(1:7,meanNoJitter, '-o', 'MarkerSize' , 5, 'MarkerFaceColor','r', 'MarkerEdgeColor', 'none', 'LineWidth',2);
    axis square;
box off;
ylabel('Potential Synapses', 'FontName', 'Arial', 'FontSize', 20);
xlabel('Cell number', 'FontName', 'Arial', 'FontSize', 20);
set(gca,'XLim',[1,7], 'FontName', 'Arial', 'FontSize', 20);

%% Plot mean of all potential synapses at different jitter radius

jitter1 = load('Jitter1um_New.mat','Pot1um');
jitter5 = load('Jitter5um_New.mat', 'Pot5um');
jitter10 = load('Jitter10_New.mat','Pot10');
jitter15 = load('Jitter15um_New.mat','Pot15um');
jitter20 = load('Jitter20um_New.mat','Pot20um');

for i = 1:22 
%     mean5(i) = sum(jitter5.Pot5um(i,:))/numel(find(jitter5.Pot5um(i,:)));
%     mean10(i) = sum(jitter10.Pot10(i,:))/numel(find(jitter10.Pot10(i,:)));
%     mean15(i) = sum(jitter15.Pot15um(i,:))/numel(find(jitter15.Pot15um(i,:)));
%     mean20(i) = sum(jitter20.Pot20um(i,:))/numel(find(jitter20.Pot20um(i,:)));

temp1Index = find(jitter1.Pot1um(i,:));
temp1Points = jitter1.Pot1um(i,temp1Index);
mean1(i) = mean(temp1Points);
stDev1(i) = std(temp1Points);
clear temp1Index;
clear temp1Points;

temp5Index = find(jitter5.Pot5um(i,:));
temp5Points = jitter5.Pot5um(i,temp5Index);
mean5(i) = mean(temp5Points);
stDev5(i) = std(temp5Points);
clear temp5Index;
clear temp5Points;

temp10Index = find(jitter10.Pot10(i,:));
temp10Points = jitter10.Pot10(i,temp10Index);
mean10(i) = mean(temp10Points);
stDev10(i) = std(temp10Points);
clear temp10Index;
clear temp10Points;

temp15Index = find(jitter15.Pot15um(i,:));
temp15Points = jitter15.Pot15um(i,temp15Index);
mean15(i) = mean(temp15Points);
stDev15(i) = std(temp15Points);
clear temp15Index;
clear temp15Points;

temp20Index = find(jitter20.Pot20um(i,:));
temp20Points = jitter20.Pot20um(i,temp20Index);
mean20(i) = mean(temp20Points);
stDev20(i) = std(temp20Points);
clear temp20Index;
clear temp20Points;


end

cellsOfInterest = find(~isnan(mean5));

% figure();
% for j = 1:numel(cellsOfInterest)
%     subplot(2,4,j);
%     histogram(jitter5.Pot5um(cellsOfInterest(j),:),'BinWidth',2); 
%     hold on;
% end
% hold off;
% 
% figure();
% for j = 1:numel(cellsOfInterest)
%     subplot(2,4,j);
%     histogram(jitter10.Pot10(cellsOfInterest(j),:),'BinWidth',2); 
%     hold on;
% end
% hold off;
% 
% figure();
% for j = 1:numel(cellsOfInterest)
%     subplot(2,4,j);
%     histogram(jitter15.Pot15um(cellsOfInterest(j),:),'BinWidth',2); 
%     hold on;
% end
% hold off;
% 
% figure();
% for j = 1:numel(cellsOfInterest)
%     subplot(2,4,j);
%     histogram(jitter20.Pot20um(cellsOfInterest(j),:),'BinWidth',2); 
%     hold on;
% end
% hold off;


figure;
hold on
plot(1:7, mean1(cellsOfInterest), '-or',  'MarkerSize' , 20, 'MarkerFaceColor','r' , 'LineWidth',2);
plot(1:7, mean5(cellsOfInterest), '-og',  'MarkerSize' , 20, 'MarkerFaceColor','g' , 'LineWidth',2);
plot(1:7, mean10(cellsOfInterest), '-ob',  'MarkerSize' , 20, 'MarkerFaceColor','b' , 'LineWidth',2);
plot(1:7, mean15(cellsOfInterest), '-om',  'MarkerSize' , 20, 'MarkerFaceColor','m' , 'LineWidth',2);
plot(1:7, mean20(cellsOfInterest), '-oc',  'MarkerSize' , 20, 'MarkerFaceColor','c' , 'LineWidth',2);
plot(1:7,meanNoJitter, '-ok', 'MarkerSize' , 20, 'MarkerFaceColor','k', 'LineWidth',2);
legend('1\mum Jitter','5\mum Jitter','10\mum Jitter','15\mum Jitter', '20\mum Jitter', 'No Jitter', 'FontName', 'Arial', 'FontSize', 20);

plot([1:7;1:7], [mean1(cellsOfInterest)-stDev1(cellsOfInterest); mean1(cellsOfInterest)+stDev1(cellsOfInterest)], 'Color','r','LineWidth',1);
plot([1:7;1:7], [mean5(cellsOfInterest)-stDev5(cellsOfInterest); mean5(cellsOfInterest)+stDev5(cellsOfInterest)], 'Color','g','LineWidth',1);
plot([1:7;1:7], [mean10(cellsOfInterest)-stDev10(cellsOfInterest); mean10(cellsOfInterest)+stDev10(cellsOfInterest)], 'Color','g','LineWidth',1);
plot([1:7;1:7], [mean15(cellsOfInterest)-stDev15(cellsOfInterest); mean15(cellsOfInterest)+stDev15(cellsOfInterest)], 'Color','b','LineWidth',1);
plot([1:7;1:7], [mean20(cellsOfInterest)-stDev20(cellsOfInterest); mean20(cellsOfInterest)+stDev20(cellsOfInterest)], 'Color','c','LineWidth',1);


axis square;
box off;
ylabel('Average potential Synapses', 'FontName', 'Arial', 'FontSize', 20);
xlabel('Cell number', 'FontName', 'Arial', 'FontSize', 20);
set(gca,'XLim',[1,7], 'FontName', 'Arial', 'FontSize', 20);
set(gcf,'color','w');

