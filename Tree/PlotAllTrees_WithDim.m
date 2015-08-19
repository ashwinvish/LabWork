
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
clc;
clear all;

% Axonal nodes for all the cells
Int1_1_axon = [];
Int1_2_axon = [];%[15 5 6 40]
Int1_3_axon = [];
Int1_4_axon = [30 42 43 44 45 47 48 41 52 49 50 51 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74];
Int1_5_axon = [56 106 107 109 110 114 115 116 159 160 161 162 163 164 165 166 167 169 170 171 168 77 80 81 82 84 86 87 85 92 99 100 103 115 117 123 125 126 137 139 150 151 153 127 124 128 131 132 134 140 143 152 154 155 156 157 158 136 138 141 142 144 106 110 114 116 133 135 145 146 147 148 149 ];
Int1_6_axon = [20 82  96 97 98 133 138 143 147 149 150 156 154 148 151 152 123  126 132 128 129 130 134 140 142 153 165 139 141 155 157 159 160 163 162 158 164 166 131 120 135 145 161 167 168 169 170 171 172 129 136 144 146 137 74 118 99];
Int1_7_axon = [9 30 75 81 99 115 120 128 110 114 116 117 121 122 123 76 107]; 
Int2_1_axon = [];%[2 5 7 14 15 25 26 73];
Int2_2_axon = [];%[5 84];
Int2_3_axon = [];%[3 67];
Int2_4_axon = [];%[2 101];
Int2_5_axon = [];%[6 64 65 66 67 69];
Int2_6_axon = [];
Int2_7_axon = [];
Int2_8_axon = [];%[10 12 22 64];
Int2_9_axon = [7 17 43 47 55 66 67 78 79 80 82 83 84 85 86 87 88 81];
Int3_1_axon = [];
Int3_2_axon = [];
Int3_3_axon = [];
Int3_4_axon = [];
Int3_5_axon = [34 52 80 82 83 86 88 90 103 104 107 111 112 113 114 115 92 81 84 87 89 95 96 98 99 100 93 94 101 102 105 106 108 109 110];
Int3_6_axon = [5 29 37 39 42 52 53 5];
MauthnerCell = [5*14474,5*49530,-45*448]; % cartesian coordinates for the center of the Mauthner cell

% all the cell IDS
cellIDs = {'Int1_1','Int1_2', 'Int1_3','Int1_4', 'Int1_5' ,'Int1_6','Int1_7' ,'Int2_1' , 'Int2_2','Int2_3' , 'Int2_4','Int2_5','Int2_6', 'Int2_7', 'Int2_8',  'Int2_9', 'Int3_1','Int3_2', 'Int3_3' 'Int3_4', 'Int3_5',  'Int3_6' };
% all the alx cells
cellIDsAlx = {'Int1_4','Int1_5','Int1_6','Int1_7','Int2_6','Int2_9','Int3_5','Int3_6'};
% all bdx1b cells
cellIDsDbx = {'Int1_2','Int1_3','Int2_1','Int2_2','Int2_3','Int2_4','Int2_5','Int2_8'};
% all barhl1 cells
cellIDsL = {'Int1_1', 'Int2_7', 'Int3_1', 'Int3_2', 'Int3_3', 'Int3_4'};

% Convert from .swc file to tree structre with presynapses,postsynapses,
% rawlegth
for kk = 1: numel(cellIDs)
    disp([cellIDs{kk} , '_WithTags.swc']);
    [thisTree,rawLength,thisPreSynapse] = generateIrreducibleDoubleLinkedTree_WithDim([cellIDs{kk} , '_WithTags.swc'],[-1:10],6, true);
    [thisTree,rawLength,thisPostSynapse] = generateIrreducibleDoubleLinkedTree_WithDim([cellIDs{kk} , '_WithTags.swc'],[-1:10],5, true);
    [thisTree,rawLength,thisSpine] = generateIrreducibleDoubleLinkedTree_WithDim([cellIDs{kk} , '_WithTags.swc'],[-1:10],9, true);
    allTrees{kk} = thisTree; allPreSynapse{kk} = thisPreSynapse; allPostSynapse{kk} = thisPostSynapse;allSpine{kk} = thisSpine;
    allRawLength{kk} = rawLength; allPost{kk} = vertcat(thisPostSynapse, thisSpine);
end

% control cells
cellControl = {'C1','C2','C3','C4','C5','C6','C7'};
TreeC1 = generateIrreducibleDoubleLinkedTree_WithDim('C1 [treeline] #317070.swc',[-1:10],5,true);
TreeC2 = generateIrreducibleDoubleLinkedTree_WithDim('C2 [treeline] #317065.swc',[-1:10],5,true);
TreeC3 = generateIrreducibleDoubleLinkedTree_WithDim('C3 [treeline] #317067.swc',[-1:10],5,true);
TreeC4 = generateIrreducibleDoubleLinkedTree_WithDim('C4 [treeline] #317072.swc',[-1:10],5,true);
TreeC5 = generateIrreducibleDoubleLinkedTree_WithDim('C5 [treeline] #317049.swc',[-1:10],5,true);
TreeC6 = generateIrreducibleDoubleLinkedTree_WithDim('C6 [treeline] #317054.swc',[-1:10],5,true);
TreeC7 = generateIrreducibleDoubleLinkedTree_WithDim('C7 [treeline] #316984.swc',[-1:10],5,true);

%% Display Control Cells

treeVisualizer(TreeC1, [1],[],[],true,{[rand rand rand]}, 1:numel(TreeC1), false);
ControlCellSoma(1,:) =  TreeC1{1,1}{1,3};
treeVisualizer(TreeC2, [1],[],[],false,{[rand rand rand]}, 1:numel(TreeC2), false);
ControlCellSoma(2,:) =  TreeC2{1,1}{1,3};
treeVisualizer(TreeC3, [1],[],[],false,{[rand rand rand]}, 1:numel(TreeC3), false);
ControlCellSoma(3,:) =  TreeC3{1,1}{1,3};
treeVisualizer(TreeC4, [1],[],[],false,{[rand rand rand]}, 1:numel(TreeC4), false);
ControlCellSoma(4,:) =  TreeC4{1,1}{1,3};

hold on;
stripe1 = [ControlCellSoma(1,:);ControlCellSoma(2,:);ControlCellSoma(3,:);ControlCellSoma(4,:)];
line(stripe1(:,1),stripe1(:,2),-stripe1(:,3),'LineWidth',2,'LineStyle','--');

treeVisualizer(TreeC5, [1],[],[],false,{[rand rand rand]}, 1:numel(TreeC5), false);
ControlCellSoma(5,:) =  TreeC5{1,1}{1,3};
treeVisualizer(TreeC6, [1],[],[],false,{[rand rand rand]}, 1:numel(TreeC6), false);
ControlCellSoma(6,:) =  TreeC6{1,1}{1,3};
treeVisualizer(TreeC7, [1],[],[],false,{[rand rand rand]}, 1:numel(TreeC7), false);
ControlCellSoma(7,:) =  TreeC7{1,1}{1,3};

 stripe2 = [ControlCellSoma(7,:);ControlCellSoma(5,:);ControlCellSoma(6,:)];
 line(stripe2(:,1),stripe2(:,2),-stripe2(:,3),'LineWidth',2,'LineStyle','--');

%%
for kk = 1:numel(cellIDs)
    % subplot(3,8,kk);
    if ismember(cellIDs{kk},cellIDsAlx)==1
        treeVisualizer(allTrees{kk}, [1],[eval([cellIDs{kk},'_axon'])],[allPost(kk) allPreSynapse(kk)],false,{[1,0.5,0]}, 1:numel(thisTree), false); % Orange for ipsiaxon
        %treeVisualizer(thisTree, [1],[],[],false,{[1,0.5,0]}, 1:numel(thisTree), false); % Orange for Alx
    elseif ismember(cellIDs{kk},cellIDsDbx)==1
        treeVisualizer(allTrees{kk}, [1],[eval([cellIDs{kk},'_axon'])],[allPost(kk) allPreSynapse(kk)],false,{[1, 0, 1]}, 1:numel(thisTree), false); % Magenta for contraxon
        %treeVisualizer(thisTree, [1],[],[],false,{[1, 0, 1]}, 1:numel(thisTree), false); % Magenta for Dbx
    else
        treeVisualizer(allTrees{kk}, [1],[eval([cellIDs{kk},'_axon'])],[allPost(kk) allPreSynapse(kk)],false,{[0, 0.5, 1]}, 1:numel(thisTree), false); % Blue for neither
        %treeVisualizer(thisTree, [1],[],[],false,{[0, 0.5, 1]}, 1:numel(thisTree), false); % Blue for lateral
    end
end
hold on;
scatter3(MauthnerCell(1,1),MauthnerCell(1,2),MauthnerCell(1,3), 50,'MarkerFaceColor','k', 'MarkerEdgeColor', 'k'); % location of the Mauther Cell center
line(stripe1(:,1),stripe1(:,2),-stripe1(:,3),'LineWidth',2,'LineStyle','--'); % stripe1
line(stripe2(:,1),stripe2(:,2),-stripe2(:,3),'LineWidth',2,'LineStyle','--'); % stripe2

h1 = gcf;
h2 = PlotViews(h1); % plots the three views of the figure handle h1
%% Plot all Soma with Time Constants
for kk = 1 : length(allTrees)
    CellSoma(kk,:) =  allTrees{1,kk}{1,1}{1,4}{1};
end

%time constants and persistence measure
tau = [7.16323454600000,6.23690945500000,9.88816223500000,17.7411357900000,5.48599643300000,10.1172859600000,7.51255049400000,18.1882559300000,57.0990589000000,100,38.7012773800000,11.6985985200000 ...
    11.0576103400000,1.44892223200000,19.5401836000000,100,16.5597699400000,10.3971051900000,26.7744003700000,7.92791676600000,6.10909614600000,5.21884052300000];
rho = [0.5207986709, 0.6661018539, 0.6648735491, 0.8190294252, 0.6370840437, 0.3520511176, 0.6180659652, 0.7770252108, 0.8445793478, 0.5048882377, 0.6169858888, 0.4693437009, 0.7319155263 ...
    0.2311146742, 0.5600555905, 0.69884915, 0.8217790613, 0.5105536878, 0.8749909893, 0.3470572669, 0.5251032841, 0.5768530374];

figure(8);
scatter3(CellSoma(:,1),CellSoma(:,2),-CellSoma(:,3), 50, rho, 'fill', 'MarkerEdgeColor', 'k');
hold on;
scatter3(MauthnerCell(1,1),MauthnerCell(1,2),MauthnerCell(1,3), 50,'MarkerFaceColor','k', 'MarkerEdgeColor', 'k');

axis([20000 140000 60000 250000]); % fixed to show the orientation of the animal
daspect([1 1 1]);
axis vis3d
box on;
set (gca,'XTick',[], 'YTick',[],'ZTick', [], 'Ydir','reverse');
view([180,90]);

figHandle = gcf;
PlotViews(figHandle);
%% ipsi,contra,unkown stats.

for kk = 1:numel(cellIDs)
    % subplot(3,8,kk);
    if ismember(cellIDs{kk},cellIDsAlx)==1
        IpsiPost(kk,1) = length(allPost{kk});
        IpsiPre(kk,1) = length(allPreSynapse{kk});
        IpsiLength(kk,1) = allRawLength{kk};
        cellClass(kk,1) = 1;
    elseif ismember(cellIDs{kk},cellIDsDbx)==1
        ContraPost(kk,1) = length(allPost{kk});
        ContraLength(kk,1) = allRawLength{kk};
        cellClass(kk,1) = 2;
    else
        UnknownPost(kk,1) = length(allPost{kk});
        UnknownLength(kk,1) = allRawLength{kk};
        cellClass(kk,1) = 3;
    end
end

% bar plot of number of Postsynapses
figure();
plot(cellClass,cellfun(@length,allPost),'o');
hold on;
meanCellPost = [mean(IpsiPost(find(IpsiPost))), mean(ContraPost(find(ContraPost))), mean(UnknownPost(find(UnknownPost))) ];
SdCellPost = [std(IpsiPost(find(IpsiPost))), std(ContraPost(find(ContraPost))), std(UnknownPost(find(UnknownPost)))];
errorbar(meanCellPost,SdCellPost,'ok');
xlabel('CellGroups');
ylabel('Number of postsynaptic sites');

figure();
plot(cellClass,cell2mat(allRawLength)/1000,'o');
hold on;
meanCellLength = [mean(IpsiLength(find(IpsiLength)))/1000, mean(ContraLength(find(ContraLength)))/1000, mean(UnknownLength(find(UnknownLength)))/1000 ]; % in um
SdCellLength = [std(IpsiLength(find(IpsiLength))/1000), std(ContraLength(find(ContraLength))/1000), std(UnknownLength(find(UnknownLength))/1000)];
errorbar(meanCellLength,SdCellLength,'ok');
xlabel('CellGroups');
ylabel('Neuronal length in \mum');


%% Get cell Somas and plot pairwise distances
% All soma locations
[Y,I] = sort(rho);
CellSomaSort = CellSoma(I,:);
tempdist = pdist(CellSoma);
EucDist = squareform(tempdist);

figure(4);
imagesc(EucDist);
colormap gray;
axis square;

Links = linkage(tempdist);
figure(5);
dendrogram(Links,'Labels',cellIDs) % dendrogram of cell distances

h = figure(6);
clear kk; clear n;
% Pairwise Eucledian distance of postsynapses
for kk = 1 : length(allTrees)
    temp  = pdist(allPost{kk});
    subplot(3,8,kk);
    imagesc(squareform(temp));
    colormap(hsv);
    axis square
    str = sprintf('%f %s', rho(kk),cellIDs{kk});
    title(str);
end
%Pairwise Eucledian distnace of Presynapse
figure(7);
clear kk; clear n;
for kk = 1 : length(allPreSynapse)
    if cellfun('isempty',allPreSynapse(1,kk)) == 1
        continue;
    else
        temp  = pdist(allPreSynapse{kk});
        subplot(3,8,kk);
        imagesc(squareform(temp));
        colormap(hsv);
        axis square
        str = sprintf('%f %s', rho(kk),cellIDs{kk});
        title(str);
    end
end

%% Additional Plots
%Rho vs Number of Post synpases
figure(11);
for i = 1:length(cellIDs)
    hold on
    if ismember(cellIDs{i},cellIDsAlx)==1
        plot(length(allPost{i}),rho(i),'Marker', '*','Color',[1 0.5 0]);
    elseif ismember(cellIDs{i},cellIDsDbx)==1
        plot(length(allPost{i}),rho(i),'Marker', '*','Color',[1 0 1]);
    else
        plot(length(allPost{i}),rho(i),'Marker', '*','Color',[0 0.5 1]);
    end
    text(length(allPost{i})+10,rho(i)+0.01,cellIDs{i});
end
%line([0,300], [0,1]);
title('Rho vs. Number of post synapses');
ylabel('Persistence measure \rho' );
xlabel('Number of postsynaptic sites');

% Rho Vs Law pathlength of neuron
figure(12)
for i = 1:length(cellIDs)
    hold on
    if ismember(cellIDs{i},cellIDsAlx)==1
        plot(cell2mat(allRawLength(1,i)),rho(i),'Marker', '*','Color',[1 0.5 0]);
    elseif ismember(cellIDs{i},cellIDsDbx)==1
        plot(cell2mat(allRawLength(1,i)),rho(i),'Marker', '*','Color',[1 0 1]);
    else
        plot(cell2mat(allRawLength(1,i)),rho(i),'Marker', '*','Color',[0 0.5 1]);
    end
    text(cell2mat(allRawLength(1,i))+1000,rho(i)+0.01,cellIDs{i});
end
%line([0,12e5],[0,1]);
title('Raw length of neuron vs. Rho');
ylabel('Persistence measure \rho');
xlabel('Raw neuron length in nm');

% Rho vs synaptic density (number of synapses/raw length)
figure(13);
for i = 1:length(cellIDs)
    hold on;
    if ismember(cellIDs{i},cellIDsAlx)==1
        plot(length(allPost{i})/cell2mat(allRawLength(1,i)),rho(i),'Marker', '*','Color',[1 0.5 0]);
    elseif ismember(cellIDs{i},cellIDsDbx)==1
        plot(length(allPost{i})/cell2mat(allRawLength(1,i)),rho(i),'Marker', '*','Color',[1 0 1]);
    else
        plot(length(allPost{i})/cell2mat(allRawLength(1,i)),rho(i),'Marker', '*','Color',[0 0.5 1]);
    end
    text(length(allPost{i})/cell2mat(allRawLength(1,i)),rho(i)+0.05,cellIDs{i});
end
title('Synaptic density of neuron vs. Rho');
ylabel('Persistence measure \rho');
xlabel('Post synaptic density (number of synapses/raw length)');

%% Gaussian Kernel
% axis([60000 250000 20000 140000 -60000 0]);
% Generate 3D gaussian Kernel

clear vol;
res = 1000; % downsampling factor
ksize = 32000; % size of Kernel in nm

% plot heat map of Postsynaptic sites
figure(14);
for kk = 1:length(cellIDs)
    subplot(3,8,kk);
    [volPost{kk},m, I] = HeatMapFish(ksize, res, allPost{kk},CellSoma(kk,:), cellIDs{kk},true);
    maxDesnity(kk) = m;
    [x,y,z] = ind2sub(size(volPost{kk}),I);
    XPost(kk) = x; YPost(kk) = y; ZPost(kk) = z-10000/res;
end

hold off;
%PostDensityPeak = [YPost'*res,XPost'*res,ZPost'*res];
PostDensityPeak = [XPost'*res,YPost'*res,ZPost'*res];

for i = 1:size(cellIDs,2)
    lengthToPostPeakNode(i) = findPathLength([cellIDs{i} , '_WithTags.swc'],[5,5,45],PostDensityPeak(i,:));
end

%% Presynapse Heatmap

figure(15);
i = 0;
for kk = 1:length(cellIDs)
    if cellfun('isempty',allPreSynapse(1,kk)) == 1
        continue;
    else
        i = i+1;
        subplot(3,8,i);
        [volPre{kk},m, I] = HeatMapFish(ksize, res, allPreSynapse{kk},CellSoma(kk,:), cellIDs{kk}, true);
    end

    maxPreDesnity(kk) = m;
    [x,y,z] = ind2sub(size(volPre{kk}),I);
    XPre(kk) = x; YPre(kk) = y; ZPre(kk) = z-10000/res;
    
end
%PreDensityPeak = [YPre'*res,XPre'*res,ZPre'*res];
PreDensityPeak = [XPre'*res,YPre'*res,ZPre'*res];

for i = 1:size(cellIDs,2)
    if cellfun('isempty',allPreSynapse(1,i)) == 1
        continue;
    else
        lengthToPrePeakNode(i) = findPathLength([cellIDs{i} , '_WithTags.swc'],[5,5,45],PreDensityPeak(i,:));
    end
end

%% dotProduct of two intersecting volumes

for i = 1:size(cellIDs,2)
    if ~cellfun('isempty',allPreSynapse(1,i)) == 1
        figure;
        for ii = 1:size(cellIDs,2)
            h = subplot(3,8,ii);
            [area] = dotVol(volPre{i},volPost{ii},CellSoma(i,:),CellSoma(ii,:),res); 
            if area == 0 
                delete (h);
            end
            IntArea(i,ii) = area;
            str = sprintf('Presynaptic cell (red): %s \nPostSynaptic cell (green): %s',cellIDs{i},cellIDs{ii});
            title(str,'FontSize',5);
        end
    else
        continue;
    end
end

%% Plot rho vs Peakdensity
%PostPeakDensity
figure(16);
for i = 1:length(cellIDs)
    hold on;
    if ismember(cellIDs{i},cellIDsAlx)==1
        plot(lengthToPostPeakNode(i),rho(i),'Marker', '*','Color',[1 0.5 0]);
    elseif ismember(cellIDs{i},cellIDsDbx)==1
        plot(lengthToPostPeakNode(i),rho(i),'Marker', '*','Color',[1 0 1]);
    else
        plot(lengthToPostPeakNode(i),rho(i),'Marker', '*','Color',[0 0.5 1]);
    end
     text(lengthToPostPeakNode(i),rho(i)+0.05,cellIDs{i});
end
xlabel('Lenght to peak post synapse density');
ylabel('rho');
%PrePeakDensity
figure(17);
for i = 1:length(cellIDs)
    hold on;
    if ismember(cellIDs{i},cellIDsAlx)==1
        plot(lengthToPrePeakNode(i),rho(i),'Marker', '*','Color',[1 0.5 0]);
    elseif ismember(cellIDs{i},cellIDsDbx)==1
        plot(lengthToPrePeakNode(i),rho(i),'Marker', '*','Color',[1 0 1]);
    else
        plot(lengthToPrePeakNode(i),rho(i),'Marker', '*','Color',[0 0.5 1]);
    end
    % text(lengthToPostPeakNode(i),rho(i)+0.05,cellIDs{i});
end
xlabel('Lenght to peak pre synapse density');
ylabel('rho')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Distribution of Synapse onto Cells

% distribution of postsynaptic sites
figure(18);
for i = 1:length(cellIDs)
    lengthToPostNode = findPathLength([cellIDs{i} , '_WithTags.swc'],[5,5,45],allPost{i});
    allLengthToPostNode{i} = lengthToPostNode ;
    subplot(3,8,i);
    scatter(1:size(allPost{i},1),sort(lengthToPostNode)/cell2mat(allRawLength(i)));
    title(cellIDs{i})
end

figure(19);
for i = 1:length(cellIDs)
    subplot(3,8,i);
    hist(sort(allLengthToPostNode{i})/cell2mat(allRawLength(i)), length(allPost{i}));title(cellIDs{i});
end

% distribution of presynaptic sites
figure(20);
for i = 1:length(cellIDs)
    if cellfun('isempty',allPreSynapse(1,i)) == 1
        continue;
    else
        lengthToPreNode = findPathLength([cellIDs{i} , '_WithTags.swc'],[5,5,45],allPreSynapse{i});
        allLengthToPreNode{i} = lengthToPreNode ;
        subplot(3,8,i);
        scatter(1:size(allPreSynapse{i},1),sort(lengthToPreNode)/cell2mat(allRawLength(i)));
        title(cellIDs{i})
    end
end
figure(21);
for i = 1:length(cellIDs)
    if cellfun('isempty',allPreSynapse(1,i)) == 1
        continue;
    else
        subplot(3,8,i);
        hist(sort(allLengthToPreNode{i})/cell2mat(allRawLength(i)), length(allPreSynapse{i}));
        title(cellIDs{i});
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Misc. plots
% plot emperical CDFs for all cells
for i = 1:22
    p = ((1:size(allLengthToPostNode{i},1))-0.5)'./size(allLengthToPostNode{i},1);
    stairs(sort(allLengthToPostNode{i})/cell2mat(allRawLength(i)),p);
    hold on;
end

% plot rawLength vs. number of PostSynapses

for i = 1:length(cellIDs)
    if ismember(cellIDs{i},cellIDsAlx)==1
        plot(cell2mat(allRawLength(i)),length(allPost{i}),'*','Color',[1 0.5 0]);
    elseif ismember(cellIDs{i},cellIDsDbx)==1
        plot(cell2mat(allRawLength(i)),length(allPost{i}),'*','Color',[1 0 1]);
    else
        plot(cell2mat(allRawLength(i)),length(allPost{i}),'*','Color',[0 0.5 1]);
    end
    hold on;
    %text(allRawLength(i),length(allPost{i}),cellIDs{i});
end

% plot cellDepth vs. Rho

for i = 1:length(cellIDs)
    if ismember(cellIDs{i},cellIDsAlx)==1
        plot(CellSoma(i,3),rho(i),'*','Color',[1 0.5 0]);
    elseif ismember(cellIDs{i},cellIDsDbx)==1
        plot(CellSoma(i,3),rho(i),'*','Color',[1 0 1]);
    else
        plot(CellSoma(i,3),rho(i),'*','Color',[0 0.5 1]);
    end
    hold on;
    %text(allRawLength(i),length(allPost{i}),cellIDs{i});
end

for i = 1:length(cellIDs)
    if ismember(cellIDs{i},cellIDsAlx)==1
        plot(i,mean(allLengthToPostNode{i}/cell2mat(allRawLength(i))) ,'*','Color',[1 0.5 0]);
        
    elseif ismember(cellIDs{i},cellIDsDbx)==1
        plot(i,mean(allLengthToPostNode{i}/cell2mat(allRawLength(i))) ,'*','Color',[1 0 1]);
    else
        plot(i,mean(allLengthToPostNode{i}/cell2mat(allRawLength(i))) ,'*','Color',[0 0.5 1]);
    end
    hold on;
    %text(allRawLength(i),length(allPost{i}),cellIDs{i});
end

%% Distribution of pre and post synaptic path lenghts
allPostSynapticLength = [];
allPreSynapticLength = [];

for i = 1:size(cellIDs,2)
allPostSynapticLength = [allPostSynapticLength;allLengthToPostNode{i}];
end

for i = 1:size(cellIDs,2)
allPreSynapticLength = [allPreSynapticLength;allLengthToPreNode{i}];
end

figure;
subplot(121);
histogram(allPostSynapticLength/1000); % dimensions in microns
title('Distribution of Postsynaptic pathlenght');
xlabel('Postsynaptic pathlenght in \mum');
ylabel('count');
subplot(122);
histogram(allPreSynapticLength/1000); % dimensions in microns
title('Distribution of Presynaptic pathlenght');
xlabel('Presynaptic pathlenght in \mum');
ylabel('count');

%% Distribution of inter-synaptic distance

clear PostSynapticDistance;
PostSynapticDistance = [];
figure;
subplot(1,2,1);
 for ii = 1:size(cellIDs,2)
     [y,I] = sort(allLengthToPostNode{ii});
     for i = 2:length(I)
         tempPost(ii,i-1) = abs(allLengthToPostNode{ii}(I(i-1)) - allLengthToPostNode{ii}(I(i)));
     end
      PostSynapticDistance = [PostSynapticDistance;tempPost(ii,:)'];
 end
 histogram(find(PostSynapticDistance>0)/1000);
 title('Inter-postsynaptic distance');
 xlabel('Distance in \mum');
 ylabel('Count');
 
 clear PreSynapticDistance;
 PreSynapticDistance = [];
 tempPre =[];
 subplot(1,2,2);
 index =1;
 for ii = 1:size(cellIDs,2)
     if ~isempty(allLengthToPreNode{ii})==1 & size(allLengthToPreNode{ii},1)>1
         [y,I] = sort(allLengthToPreNode{ii});
         for i = 2:length(I)
             tempPre(ii,i-1) = abs(allLengthToPreNode{ii}(I(i-1)) - allLengthToPreNode{ii}(I(i)));
         end
         PreSynapticDistance = [PreSynapticDistance;tempPre(index,:)'];
         index = index+1;
     else
         index = index+1;
         continue;
     end
 end
 histogram(find(PreSynapticDistance>0)/1000);
 title('Inter-presynaptic distance');
 xlabel('Distance in \mum');
 ylabel('Count');
