for i = 1:size(cellIDs,2)
    treeVisualizer(allTrees{i}, [1],[eval([cellIDs{i},'_axon'])],[{allPostSynapse{i}} {allPreSynapse{i}}],false,{[rand,rand,rand]}, 1:numel(allTrees{i}), false);
    h1 = gcf;
    h2 = PlotViews(h1);
    pause(10);
    saveas(gcf,cellIDs{i},'png');
    close all;
end



%%


ksize = 1000*[2:6:52];
%area = zeros(size(ksize));
postCell = [2,9,21];
clear area;

for kk = 1:length(postCell)
    figure(kk);
    ha = tight_subplot(2,5,[.05 .05],[.05 .1],[.01 .01]);
    for i = 1:length(ksize)
        vol1 = HeatMapFish(ksize(i), 1000, allPreSynapse{5},CellSoma(5,:), cellIDs{5},false);
        vol2 = HeatMapFish(ksize(i), 1000, allPostSynapse{postCell(kk)},CellSoma(clc
        linpostCell(kk),:), cellIDs{postCell(kk)},false);
        axes(ha(i));
        area(kk,i) = dotVol(vol1,vol2,CellSoma(5,:),CellSoma(postCell(kk),:),1000);
        str1 = sprintf('Kernel Size %d',ksize(i)/1000);
        title(str1);
        str2 = sprintf('Presynaptic cell (red): %s \nPostSynaptic cell (green): %s',cellIDs{5},cellIDs{postCell(kk)});
        figtitle(str2);
        hold on;
    end
    delete(ha(10));
end

figure;

plot(ksize/1000,area(1,:)/10e6,'-*');
hold on;
plot(ksize/1000,area(2,:)/10e6,'-*');
plot(ksize/1000,area(3,:)/10e6,'-*');
xlabel('Kernel size in \mum');
ylabel('Overlap area in \mum^2');

%%
ha = tight_subplot(2,5,[.05 .05],[.05 .1],[.01 .01]);
for i = 1:size(ksize,2)
    vol1 = HeatMapFish(ksize(i), 1000, allPreSynapse{5},CellSoma(5,:), cellIDs{5},false);
    vol2 = HeatMapFish(ksize(i), 1000, allPostSynapse{2},CellSoma(2,:), cellIDs{2},false);
    axes(ha(i));
    area(kk,i) = dotVol(vol1,vol2,CellSoma(5,:),CellSoma(postCell(kk),:),1000);
    str1 = sprintf('Kernel Size %d',ksize(i)/1000);
        title(str1);
        str2 = sprintf('Presynaptic cell (red): %s \nPostSynaptic cell (green): %s',cellIDs{5},cellIDs{postCell(kk)});
        figtitle(str2);
end

%% ratio of number of post/pre synapses
numberOfPresynapses = [];
numberOfPostsynapses = [];
for i = 1:size(cellIDs,2)
    numberOfPresynapses = [ numberOfPresynapses, length(allPreSynapse{i})];
    numberOfPostsynapses = [ numberOfPostsynapses, length(allPostSynapse{i})];
end

%sprintf('Number of postsynapses/ number of presynapses = %3d', numberOfPostsynapses./numberOfPresynapses)

%% ratio of dendritic length/ axonal length
axLength = [];
clear temp;
for i = 1:size(cellIDs,2)
    if eval([cellIDs{i},'_axon'])>0
           temp(1:size(allLengthToPreNode{i},1)-1,i) = diff(sort(allLengthToPreNode{i}));
           axLength = [axLength,sum(temp(:,i))];
    else
        axLength = [axLength,0];
        continue;
    end
end
denLength = cell2mat(allRawLength)-axLength;
sprintf('dendrite length / axon length = %d',sum(denLength)/sum(axLength))

%% ratio of dendritic length/ axonal length v2

%% Euclidean distance from Mauthner Cell


dist = sqrt((MauthnerCell(1,1)-CellSoma(:,1)).^2 + (MauthnerCell(1,2)-CellSoma(:,2)).^2 + (MauthnerCell(1,3)-CellSoma(:,3)).^2);
dist = dist/1000 ; % in um
distRC = sqrt((MauthnerCell(1,1)-CellSoma(:,1)).^2 + (MauthnerCell(1,2)-CellSoma(:,2)).^2);
distRC = distRC/1000;% in um
distDV  =  sqrt((MauthnerCell(1,1)-CellSoma(:,1)).^2);
distDV = distDV/1000;% in um

figure(1)
for i = 1:length(cellIDs)
    if ismember(cellIDs{i},cellIDsAlx)==1
        plot(dist(i),rho(i),'Marker', '*','Color',[1 0.5 0]);
        hold on;
    elseif ismember(cellIDs{i},cellIDsDbx)==1
        plot(dist(i),rho(i),'Marker', '*','Color',[1 0 1]);
    else
        plot(dist(i),rho(i),'Marker', '*','Color',[0 0.5 1]);;
    end
end
title('Distance from center of Mauthner cell');
xlabel('Euclidean distance to Mauthner cell in \mum');
ylabel('Peristance measuere \rho');

figure(2);
for i = 1:length(cellIDs)
    if ismember(cellIDs{i},cellIDsAlx)==1
        plot(distRC(i),rho(i),'Marker', '*','Color',[1 0.5 0]);
        hold on;
    elseif ismember(cellIDs{i},cellIDsDbx)==1    
        plot(distRC(i),rho(i),'Marker', '*','Color',[1 0 1]);     
    else
        plot(distRC(i),rho(i),'Marker', '*','Color',[0 0.5 1]);
    end
end
title('Rostro-Caudal distance from center of Mauthner cell');
xlabel('Rostro-Caudal Euclidean distance to Mauthner cell in \mum');
ylabel('Peristance measuere \rho');

figure(3);
for i = 1:length(cellIDs)
    if ismember(cellIDs{i},cellIDsAlx)==1
        plot(distDV(i),rho(i),'Marker', '*','Color',[1 0.5 0]);
        hold on;
    elseif ismember(cellIDs{i},cellIDsDbx)==1 
        plot(distDV(i),rho(i),'Marker', '*','Color',[1 0 1]);
    else   
        plot(distDV(i),rho(i),'Marker', '*','Color',[0 0.5 1]);
    end
end

title('Dorso-ventral distance from center of Mauthner cell');
xlabel('Dorso-ventral Euclidean distance to Mauthner cell in \mum');
ylabel('Peristance measuere \rho');

    