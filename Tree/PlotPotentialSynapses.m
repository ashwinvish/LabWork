% %%
PotSites = [];
PotSitesTemp = [];
index = 1;
cellIDAlxSparse = {'Int1_4','Int1_5','Int1_6','Int1_7','Int2_6','Int2_9','Int3_5','Int3_6'};

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
plot(i,s,'o', 'MarkerSize' , 25, 'MarkerFaceColor','k', 'MarkerEdgeColor', 'none');
hold on;
meanNoJitter = sumOfPotentialSynapses./nonZeroElements;
meanNoJitter(1) = 0;
plot(1:8,meanNoJitter, 'o', 'MarkerSize' , 35, 'MarkerFaceColor','r', 'MarkerEdgeColor', 'none', 'LineWidth',2);
axis square;
box off;
ylabel('Potential synapses', 'FontName', 'Arial', 'FontSize', 40);
xlabel('Neuron #', 'FontName', 'Arial', 'FontSize', 40);
set(gca,'XLim',[1,8], 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
set(gcf,'color','w');

%% Plot mean of all potential synapses at different jitter radius

jitter1 = load('Shuffle1000um_06-Feb-2016.mat','PotSites');
jitter5 = load('Shuffle5000um_20-Jan-2016.mat', 'PotSites');
jitter10 = load('Shuffle10000um_20-Jan-2016.mat','PotSites');
jitter15 = load('Shuffle15000um_21-Jan-2016.mat','PotSites');
jitter20 = load('Shuffle20000um_21-Jan-2016.mat','PotSites');

for i = 1:22
    %     mean5(i) = sum(jitter5.Pot5um(i,:))/numel(find(jitter5.Pot5um(i,:)));
    %     mean10(i) = sum(jitter10.Pot10(i,:))/numel(find(jitter10.Pot10(i,:)));
    %     mean15(i) = sum(jitter15.Pot15um(i,:))/numel(find(jitter15.Pot15um(i,:)));
    %     mean20(i) = sum(jitter20.Pot20um(i,:))/numel(find(jitter20.Pot20um(i,:)));
    
    temp1Index = find(jitter1.PotSites(i,:));
    temp1Points = jitter1.PotSites(i,temp1Index);
    mean1(i) = mean(temp1Points);
    stDev1(i) = std(temp1Points);
    clear temp1Index;
    clear temp1Points;
    
    temp5Index = find(jitter5.PotSites(i,:));
    temp5Points = jitter5.PotSites(i,temp5Index);
    mean5(i) = mean(temp5Points);
    stDev5(i) = std(temp5Points);
    clear temp5Index;
    clear temp5Points;
    
    temp10Index = find(jitter10.PotSites(i,:));
    temp10Points = jitter10.PotSites(i,temp10Index);
    mean10(i) = mean(temp10Points);
    stDev10(i) = std(temp10Points);
    clear temp10Index;
    clear temp10Points;
    
    temp15Index = find(jitter15.PotSites(i,:));
    temp15Points = jitter15.PotSites(i,temp15Index);
    mean15(i) = mean(temp15Points);
    stDev15(i) = std(temp15Points);
    clear temp15Index;
    clear temp15Points;
    
    temp20Index = find(jitter20.PotSites(i,:));
    temp20Points = jitter20.PotSites(i,temp20Index);
    mean20(i) = mean(temp20Points);
    stDev20(i) = std(temp20Points);
    clear temp20Index;
    clear temp20Points;
    
    
end

cellsOfInterest = [4,5,6,7,13,16,21,22]; % both Alx and Trans cells


figure();
for j = 1:numel(cellsOfInterest)
    subplot(3,4,j);
    histogram(jitter1.PotSites(cellsOfInterest(j),:),'BinWidth',1);
    set(gca,'XLim',[0, max(jitter1.PotSites(cellsOfInterest(j),:))+1]);
    title(cellIDs(cellsOfInterest(j)));
    hold on;
end
subplot(3,4,9);
plot(1:8, mean1(cellsOfInterest), '-or',  'MarkerSize' , 10, 'MarkerFaceColor','r' , 'LineWidth',2);
hold on;
plot(1:8,sumOfPotentialSynapses, '-ok', 'MarkerSize' , 10, 'MarkerFaceColor','k', 'LineWidth',2);
plot([1:8;1:8], [mean1(cellsOfInterest)-stDev1(cellsOfInterest); mean1(cellsOfInterest)+stDev1(cellsOfInterest)], 'Color','r','LineWidth',1);
ylabel('Average potential Synapses');
xlabel('Cell number');
set(gcf,'color','w');
suptitle('JitterRadius 1\mum');
hold off;

figure();
for j = 1:numel(cellsOfInterest)
    subplot(3,4,j);
    histogram(jitter5.PotSites(cellsOfInterest(j),:),'BinWidth',1);
    set(gca,'XLim',[0, max(jitter5.PotSites(cellsOfInterest(j),:))+1]);
    title(cellIDs(cellsOfInterest(j)));
    hold on;
end
subplot(3,4,9);
plot(1:8, mean5(cellsOfInterest), '-og',  'MarkerSize' , 10, 'MarkerFaceColor','g' , 'LineWidth',2);
hold on;
plot(1:8,sumOfPotentialSynapses, '-ok', 'MarkerSize' , 10, 'MarkerFaceColor','k', 'LineWidth',2);
plot([1:8;1:8], [mean5(cellsOfInterest)-stDev5(cellsOfInterest); mean5(cellsOfInterest)+stDev5(cellsOfInterest)], 'Color','g','LineWidth',1);
ylabel('Average potential Synapses');
xlabel('Cell number');
set(gcf,'color','w');
suptitle('JitterRadius 5\mum');
hold off;

figure();
for j = 1:numel(cellsOfInterest)
    subplot(3,4,j);
    histogram(jitter10.PotSites(cellsOfInterest(j),:),'BinWidth',1);
    set(gca,'XLim',[0, max(jitter10.PotSites(cellsOfInterest(j),:))+1]);
    title(cellIDs(cellsOfInterest(j)));
    hold on;
end
subplot(3,4,9);
plot(1:8, mean10(cellsOfInterest), '-ob',  'MarkerSize' , 10, 'MarkerFaceColor','b' , 'LineWidth',2);
hold on;
plot(1:8,sumOfPotentialSynapses, '-ok', 'MarkerSize' , 10, 'MarkerFaceColor','k', 'LineWidth',2);
plot([1:8;1:8], [mean10(cellsOfInterest)-stDev10(cellsOfInterest); mean10(cellsOfInterest)+stDev10(cellsOfInterest)], 'Color','b','LineWidth',1);
ylabel('Average potential Synapses');
xlabel('Cell number');
set(gcf,'color','w');
suptitle('JitterRadius 10\mum');
hold off;

figure();
for j = 1:numel(cellsOfInterest)
    subplot(3,4,j);
    histogram(jitter15.PotSites(cellsOfInterest(j),:),'BinWidth',1);
    set(gca,'XLim',[0, max(jitter15.PotSites(cellsOfInterest(j),:))+1]);
    title(cellIDs(cellsOfInterest(j)));
    hold on;
end
subplot(3,4,9);
plot(1:8, mean15(cellsOfInterest), '-om',  'MarkerSize' , 10, 'MarkerFaceColor','m' , 'LineWidth',2);
hold on;
plot(1:8,sumOfPotentialSynapses, '-ok', 'MarkerSize' , 10, 'MarkerFaceColor','k', 'LineWidth',2);
plot([1:8;1:8], [mean15(cellsOfInterest)-stDev15(cellsOfInterest); mean15(cellsOfInterest)+stDev15(cellsOfInterest)], 'Color','m','LineWidth',1);
ylabel('Average potential Synapses');
xlabel('Cell number');
set(gcf,'color','w');
suptitle('JitterRadius 15\mum');
hold off;

figure();
for j = 1:numel(cellsOfInterest)
    subplot(3,4,j);
    histogram(jitter20.PotSites(cellsOfInterest(j),:),'BinWidth',1);
    set(gca,'XLim',[0, max(jitter20.PotSites(cellsOfInterest(j),:))+1]);
    title(cellIDs(cellsOfInterest(j)));
    hold on;
end
subplot(3,4,9);
plot(1:8, mean20(cellsOfInterest), '-oc',  'MarkerSize' , 10, 'MarkerFaceColor','c' , 'LineWidth',2);
hold on;
plot(1:8,sumOfPotentialSynapses, '-ok', 'MarkerSize' , 10, 'MarkerFaceColor','k', 'LineWidth',2);
plot([1:8;1:8], [mean20(cellsOfInterest)-stDev20(cellsOfInterest); mean20(cellsOfInterest)+stDev20(cellsOfInterest)], 'Color','c','LineWidth',1);
ylabel('Average potential Synapses');
xlabel('Cell number');
set(gcf,'color','w');
suptitle('JitterRadius 20\mum');
hold off;


figure;
hold on
plot(1:8, mean1(cellsOfInterest), 'or',  'MarkerSize' , 35, 'MarkerFaceColor','r' , 'LineWidth',2);
plot(1:8, mean5(cellsOfInterest), 'og',  'MarkerSize' , 35, 'MarkerFaceColor','g' , 'LineWidth',2);
plot(1:8, mean10(cellsOfInterest), 'ob',  'MarkerSize' , 35, 'MarkerFaceColor','b' , 'LineWidth',2);
plot(1:8, mean15(cellsOfInterest), 'om',  'MarkerSize' , 35, 'MarkerFaceColor','m' , 'LineWidth',2);
plot(1:8, mean20(cellsOfInterest), 'oc',  'MarkerSize' , 35, 'MarkerFaceColor','c' , 'LineWidth',2);
plot(1:8, sumOfPotentialSynapses, 'ok', 'MarkerSize' , 35, 'MarkerFaceColor','k', 'LineWidth',2);

legend({'1\mum ','5\mum ','10\mum ','15\mum ', '20\mum ', 'No Jitter'},'FontName', 'Arial', 'FontSize', 40,'Box','off');

plot([1:8;1:8], [mean1(cellsOfInterest)-stDev1(cellsOfInterest); mean1(cellsOfInterest)+stDev1(cellsOfInterest)], 'Color','r','LineWidth',2);
plot([1:8;1:8], [mean5(cellsOfInterest)-stDev5(cellsOfInterest); mean5(cellsOfInterest)+stDev5(cellsOfInterest)], 'Color','g','LineWidth',2);
plot([1:8;1:8], [mean10(cellsOfInterest)-stDev10(cellsOfInterest); mean10(cellsOfInterest)+stDev10(cellsOfInterest)], 'Color','b','LineWidth',2);
plot([1:8;1:8], [mean15(cellsOfInterest)-stDev15(cellsOfInterest); mean15(cellsOfInterest)+stDev15(cellsOfInterest)], 'Color','m','LineWidth',2);
plot([1:8;1:8], [mean20(cellsOfInterest)-stDev20(cellsOfInterest); mean20(cellsOfInterest)+stDev20(cellsOfInterest)], 'Color','c','LineWidth',2);


axis square;
box off;
ylabel('Potential synapses', 'FontName', 'Arial', 'FontSize', 40);
xlabel('Neuron #', 'FontName', 'Arial', 'FontSize', 40);
set(gca,'XLim',[1,8], 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
set(gcf,'color','w');
%set(gca,'color','none');

%%
figure()
hold on;
JitterSteps = [1,5,10,15,20];
mean0 = sumOfPotentialSynapses;
for i = 1:numel(cellsOfInterest)
    for j = 1:numel(JitterSteps)
        str = sprintf('mean%d',JitterSteps(j));
        plot(j,str(cellsOfInterest(i)),'o');
        hold on;
    end
end
% plot(1, mean1(cellsOfInterest), '-or',  'MarkerSize' , 20, 'MarkerFaceColor','r' , 'LineWidth',2);
% plot(2, mean5(cellsOfInterest), '-og',  'MarkerSize' , 20, 'MarkerFaceColor','g' , 'LineWidth',2);
% plot(3, mean10(cellsOfInterest), '-ob',  'MarkerSize' , 20, 'MarkerFaceColor','b' , 'LineWidth',2);
% plot(4, mean15(cellsOfInterest), '-om',  'MarkerSize' , 20, 'MarkerFaceColor','m' , 'LineWidth',2);
% plot(5, mean20(cellsOfInterest), '-oc',  'MarkerSize' , 20, 'MarkerFaceColor','c' , 'LineWidth',2);
% plot(0, sumOfPotentialSynapses, '-ok',  'MarkerSize' , 20, 'MarkerFaceColor','k' , 'LineWidth',2);

axis square;
box off;
ylabel('Average potential Synapses', 'FontName', 'Arial', 'FontSize', 20);
xlabel('Jitter', 'FontName', 'Arial', 'FontSize', 20);
set(gca,'XLim',[0,5], 'FontName', 'Arial', 'FontSize', 20);
set(gcf,'color','w');

%%

