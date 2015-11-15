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
            if isempty(validNodes)                                                   % if no axon check
                validNodes = 1:numel(DenTree);
            end                                                                      % if no axon then iterate through all nodes
            for jj = 1:numel(DenTree)
                children = DenTree{jj}{2};                                           % Consider childern of jj node
                for nn = 1:numel(children)                                           % iterate over all children
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
    %suptitle(cellIDs(i));
    
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
    plot(1:7,meanNoJitter, '-o', 'MarkerSize' , 5, 'MarkerFaceColor','r', 'MarkerEdgeColor', 'none');
    axis square;
box off;
ylabel('Potential Synapses', 'FontName', 'Arial', 'FontSize', 20);
xlabel('Cell number', 'FontName', 'Arial', 'FontSize', 20);
set(gca,'XLim',[1,7], 'FontName', 'Arial', 'FontSize', 20);

%%



jitter5 = load('Jitter5um_New.mat', 'Pot5um');
jitter10 = load('Jitter10_New.mat','Pot10');
jitter15 = load('Jitter15um_New.mat','Pot15um');
jitter20 = load('Jitter20um_New.mat','Pot20um');

for i = 1:22 
    mean5(i) = sum(jitter5.Pot5um(i,:))/numel(find(jitter5.Pot5um(i,:)));
    mean10(i) = sum(jitter10.Pot10(i,:))/numel(find(jitter10.Pot10(i,:)));
    mean15(i) = sum(jitter15.Pot15um(i,:))/numel(find(jitter15.Pot15um(i,:)));
    mean20(i) = sum(jitter20.Pot20um(i,:))/numel(find(jitter20.Pot20um(i,:)));
end

cellsOfInterest = find(~isnan(mean5));

figure();
subplot(2,2,1);
for j = 1:numel(cellsOfInterest)
    histogram(jitter5.Pot5um(cellsOfInterest(j),:),'BinWidth',2); 
    hold on;
end
hold off;
subplot(2,2,2);
for j = 1:numel(cellsOfInterest)
    histogram(jitter10.Pot10(cellsOfInterest(j),:),'BinWidth',2); 
    hold on;
end
hold off;
subplot(2,2,3);
for j = 1:numel(cellsOfInterest)
    histogram(jitter15.Pot15um(cellsOfInterest(j),:),'BinWidth',2); 
    hold on;
end
hold off;
subplot(2,2,4);
for j = 1:numel(cellsOfInterest)
    histogram(jitter20.Pot20um(cellsOfInterest(j),:),'BinWidth',2); 
    hold on;
end
hold off;


figure;
hold on
plot(1:7, mean5(cellsOfInterest), '-or',  'MarkerSize' , 20, 'MarkerFaceColor','r' , 'LineWidth',2);
plot(1:7, mean10(cellsOfInterest), '-og',  'MarkerSize' , 20, 'MarkerFaceColor','g' , 'LineWidth',2);
plot(1:7, mean15(cellsOfInterest), '-ob',  'MarkerSize' , 20, 'MarkerFaceColor','b' , 'LineWidth',2);
plot(1:7, mean20(cellsOfInterest), '-oc',  'MarkerSize' , 20, 'MarkerFaceColor','c' , 'LineWidth',2);
plot(1:7,meanNoJitter, '-ok', 'MarkerSize' , 20, 'MarkerFaceColor','k', 'LineWidth',2);


axis square;
box off;
ylabel('Average potential Synapses', 'FontName', 'Arial', 'FontSize', 20);
xlabel('Cell number', 'FontName', 'Arial', 'FontSize', 20);
set(gca,'XLim',[1,8], 'FontName', 'Arial', 'FontSize', 20);
set(gcf,'color','w');

legend('5\mum Jitter','10\mum Jitter','15\mum Jitter', '20\mum Jitter', 'No Jitter', 'FontName', 'Arial', 'FontSize', 20)
%set(gca,'YLim',[0 10);


