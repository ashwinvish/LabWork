%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inducing nodes
%Int1_1 = none;1
%Int1_2 = [15 5 6 40];
%Int1_3 = [1 71];
%Int1_4 = [39 46 47 48 52 49 50 51 53 54 55 56 57 58 59 60 62 63 64 65 66 61]; 
%Int1_5 = [56 75 91 93 104 105 106 107 108 109 110 111 112 114 115 116 113 ];
%Int1_6 = [11 17 20 26 82 95 96 97 98 133 138 143 147 149 150 156 154 148 151 152 122 125 127 132 128 130 134 140 142 153 165 139 141 155 157 159 160 163 162 158 164 166 131 117 119 120 135 145 161 167 168 169 170 171 172 129 136 144 146 137 56 66 68 72 75 83 117 86 99];
%Int1_7 = [9 30 70 71 76 92 105 109 113 103 104 106 107 110 111 112];
%Int2_1 = [2 5 7 14 15 25 26 73];
%Int2_2 = [5 84];
%Int2_3 = [3 67];
%Int2_4 = [2 101];
%Int2_5 = [6 64 65 66 67 69];
%Int2_8 = [10 12 22 64];
%Int3_5 = [1 91 34 52 80 82 83 86 88 90 103 104 107 111 112 113 114 115 92 81 84 87 89 95 96 98 99 100 93 94 101 102]; 

clc; 
clear all;
% all the cell IDS
cellIDs = {'Int1_1','Int1_2', 'Int1_3','Int1_4', 'Int1_5' ,'Int1_6','Int1_7' ,'Int2_1' , 'Int2_2','Int2_3' , 'Int2_4','Int2_5','Int2_6', 'Int2_7', 'Int2_8',  'Int2_9', 'Int3_1','Int3_2', 'Int3_3' 'Int3_4', 'Int3_5',  'Int3_6' };
% all the alx cells
cellIDsAlx = {'Int1_4','Int1_5','Int1_6','Int1_7','Int2_6','Int2_9','Int3_5','Int3_6'};
% all bdx1b cells
cellIDsDbx = {'Int1_2','Int1_3','Int2_1','Int2_2','Int2_3','Int2_4','Int2_5','Int2_8','Int3_5'};
% all barhl1 cells
cellIDsL = {'Int1_1', 'Int2_7', 'Int3_1', 'Int3_2', 'Int3_3', 'Int3_4'}; 


for kk = 1: numel(cellIDs)
    disp([cellIDs{kk} , '_WithTags.swc']);
    [thisTree,rawLength,thisPreSynapse] = generateIrreducibleDoubleLinkedTree_WithDim([cellIDs{kk} , '_WithTags.swc'],[-1:10],6, true);
    [thisTree,rawLength,thisPostSynapse] = generateIrreducibleDoubleLinkedTree_WithDim([cellIDs{kk} , '_WithTags.swc'],[-1:10],5, true);
    [thisTree,rawLength,thisSpine] = generateIrreducibleDoubleLinkedTree_WithDim([cellIDs{kk} , '_WithTags.swc'],[-1:10],9, true);
    allTrees{kk} = thisTree; allPreSynapse{kk} = thisPreSynapse; allPostSynapse{kk} = thisPostSynapse;allSpine{kk} = thisSpine;
    allRawLength(kk) = rawLength; allPost{kk} = vertcat(thisPostSynapse, thisSpine);
   % subplot(3,8,kk);
    if ismember(cellIDs{kk},cellIDsAlx)==1
        treeVisualizer(thisTree, [1],[],[{thisPreSynapse} {thisPostSynapse}],false,{[1,0.5,0]}, 1:numel(thisTree), false); % Orange for Alx
        %treeVisualizer(thisTree, [1],[],[],false,{[1,0.5,0]}, 1:numel(thisTree), false); % Orange for Alx
    elseif ismember(cellIDs{kk},cellIDsDbx)==1
        %treeVisualizer(thisTree, [1],[],[{thisPreSynapse} {thisPostSynapse}],false,{[1, 0, 1]}, 1:numel(thisTree), false); % Magenta for Dbx
        %treeVisualizer(thisTree, [1],[],[],false,{[1, 0, 1]}, 1:numel(thisTree), false); % Magenta for Dbx
    else
        %treeVisualizer(thisTree, [1],[],[{thisPreSynapse} {thisPostSynapse}],false,{[0, 0.5, 1]}, 1:numel(thisTree), false); % Blue for lateral
        %treeVisualizer(thisTree, [1],[],[],false,{[0, 0.5, 1]}, 1:numel(thisTree), false); % Blue for lateral
    end
    daspect([1 1 1]);
    
    %view([-50,40]);
    % box on;
    %grid on;
    %daspect([1 1 1]);
    % axis vis3d;
    %XColor = [1,1,1]; YColor = [1,1,1];
    %h = gcf;
    %set (gca, 'XTick',[], 'YTick',[],'ZTick', []);
    %view([-90,90]);
    %title(cellIDs{kk});
    %saveas(h,[cellIDs{kk},'.png'],'png');
end
% daspect([1 1 1]);
h1 = gcf;
box on;
XColor = [1,1,1]; YColor = [1,1,1];
set (gca, 'XTick',[], 'YTick',[],'ZTick', []);
%view([-50,40]);
view([-90,90]);
%XZ view
h2 = figure(2);
copyobj(get(h1,'children'),h2);
axis([-0 30000 20000 140000]);
view([-90, 0]); 
daspect([2 1 2]);
% YZ view
h3 = figure(3);
copyobj(get(h1,'children'),h3);
axis([60000 250000 -30000 0]);
view([0,0]);camroll(90); 
%% Get cell Somas and plot pairwise distances
% All soma locations
for kk = 1 : length(allTrees)
    thisCellSoma = allTrees{1,kk}{1,1}{1,4}{1};
    cellSoma{kk} = thisCellSoma;
end
temp = cell2mat(cellSoma);
CellSoma = vec2mat(temp,3);
CellSoma(:,[1,2,3]) = CellSoma(:,[2,1,3]);
[Y,I] = sort(rho);
CellSomaSort = CellSoma(I,:);
tempdist = pdist(CellSomaSort);
EucDist = squareform(tempdist);
figure(4);
imagesc(EucDist);
colormap(hsv);
axis square;

Links = linkage(tempdist);
figure(5); 
dendrogram(Links,'Labels',cellIDs) % dendrogram of cell distances

h = figure(6);
clear kk; clear n;
% Pairwise Eucledian distance of postsynapses
for kk = 1 : length(allTrees)
    for n = 1:size(allPostSynapse{kk},1)
        X = allPostSynapse{kk}(n,1);
        Y = allPostSynapse{kk}(n,2);
        Z = allPostSynapse{kk}(n,3);
        PostSynapse(n,:) = horzcat(Y,X,Z);
    end
    allPostSynapseSwapped{kk} = PostSynapse;
    clear PostSynapse;
    temp  = pdist(allPostSynapseSwapped{kk});
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
    for n = 1:size(allPreSynapse{kk},1)
        X = allPreSynapse{kk}(n,1);
        Y = allPreSynapse{kk}(n,2);
        Z = allPreSynapse{kk}(n,3);
        PreSynapse(n,:) = horzcat(Y,X,Z);
    end
    allPreSynapseSwapped{kk} = PreSynapse;
    clear PreSynapse;
    temp  = pdist(allPreSynapseSwapped{kk});
    subplot(3,8,kk);
    imagesc(squareform(temp));
    colormap(hsv);
    axis square
    str = sprintf('%f %s', rho(kk),cellIDs{kk});
    title(str);
    end
end



%% Plot all Soma with Time Constants
tau = [7.16323454600000,6.23690945500000,9.88816223500000,17.7411357900000,5.48599643300000,10.1172859600000,7.51255049400000,18.1882559300000,57.0990589000000,100,38.7012773800000,11.6985985200000,11.0576103400000,1.44892223200000,19.5401836000000,100,16.5597699400000,10.3971051900000,26.7744003700000,7.92791676600000,6.10909614600000,5.21884052300000];
rho = [0.5207986709, 0.6661018539, 0.6648735491, 0.8190294252, 0.6370840437, 0.3520511176, 0.6180659652, 0.7770252108, 0.8445793478, 0.5048882377, 0.6169858888, 0.4693437009, 0.7319155263 ...
0.2311146742, 0.5600555905, 0.69884915, 0.8217790613, 0.5105536878, 0.8749909893, 0.3470572669, 0.5251032841, 0.5768530374];
figure(8);
scatter3(CellSoma(:,1),CellSoma(:,2),-CellSoma(:,3), 50, rho, 'fill', 'MarkerEdgeColor', 'k');
grid off;
daspect([1 1 1]);
axis([60000 250000 20000 140000]);
box on;
view([-90,90]);
XColor = [1,1,1]; YColor = [1,1,1];
hold on;
plot([240000, 240000], [120000, 140000], '-k' );
% to plot the other views
h1 = gcf;
set (gca, 'XTick',[], 'YTick',[],'ZTick', []);
h2 = figure(9);
copyobj(get(h1,'children'),h2);
view([-90, 0]); % xz view
h3 = figure(10);
copyobj(get(h1,'children'),h3);
view([0,0]);camroll(90); % yz view
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
ylabel('Persistence measure Rho' );
xlabel('Number of postsynaptic sites');

% Rho Vs Law pathlength of neuron
figure(12)
for i = 1:length(cellIDs)
    hold on
    if ismember(cellIDs{i},cellIDsAlx)==1
        plot(allRawLength(1,i),rho(i),'Marker', '*','Color',[1 0.5 0]);
    elseif ismember(cellIDs{i},cellIDsDbx)==1
        plot(allRawLength(1,i),rho(i),'Marker', '*','Color',[1 0 1]);
    else
        plot(allRawLength(1,i),rho(i),'Marker', '*','Color',[0 0.5 1]);
    end
    text(allRawLength(1,i)+1000,rho(i)+0.01,cellIDs{i});
end
%line([0,12e5],[0,1]);
title('Raw length of neuron vs. Rho');
ylabel('Persistence measure Rho');
xlabel('Raw neuron length in nm');

% Rho vs synaptic density (number of synapses/raw length)
figure(13);
for i = 1:length(cellIDs)
    hold on;
    if ismember(cellIDs{i},cellIDsAlx)==1
        plot(length(allPost{i})/allRawLength(1,i),rho(i),'Marker', '*','Color',[1 0.5 0]);
    elseif ismember(cellIDs{i},cellIDsDbx)==1
        plot(length(allPost{i})/allRawLength(1,i),rho(i),'Marker', '*','Color',[1 0 1]);
    else
        plot(length(allPost{i})/allRawLength(1,i),rho(i),'Marker', '*','Color',[0 0.5 1]);
    end
    text(length(allPost{i})/allRawLength(1,i),rho(i)+0.05,cellIDs{i});
end
title('Synaptic density of neuron vs. Rho');
ylabel('Persistence measure Rho');
xlabel('Post synaptic density ( number of synapses/raw length)');

%% Gaussian Kernel
% axis([60000 250000 20000 140000 -60000 0]);
% Generate 3D gaussian Kernel

clear vol;
res =1000; % resolution of image
ksize = 18; % size of Kernel
Gwin = gausswin(ksize+1);
K = Gwin*Gwin';

for i = 1:(ksize+1)
    K3D(:,:,i) = K(:,:)*Gwin(i);
end

% plot heat map of presynaptic sites
figure(14);
for kk = 1:length(allPost)
    volPost = zeros((15e4-0)/res,(2.5e5-0)/res,(8e4-0)/res); 
    % add 10000/res pixels to the zaxis to avoid edge effects
    for n = 1:length(allPost{kk});
        volPost(round(allPost{kk}(n,1)/res) + (-ksize/2:ksize/2), round(allPost{kk}(n,2)/res) + (-ksize/2:ksize/2), round(allPost{kk}(n,3)/res + 10000/res)+ (-ksize/2:ksize/2)) = ...
            volPost(round(allPost{kk}(n,1)/res) + (-ksize/2:ksize/2), round(allPost{kk}(n,2)/res) + (-ksize/2:ksize/2),round(allPost{kk}(n,3)/res + 10000/res)+ (-ksize/2:ksize/2))  + K3D;
    end
    subplot(3,8,kk);
    %plot XZ
    tempXZ = max(volPost,[],2);
    for i = 1:size(tempXZ,3)
        B(:,i) = tempXZ(:,:,i);
    end
    imagesc(0,0,imrotate(B,-90)); colormap(jet);
    title(cellIDs{kk});
    hold on;
    plot(160-CellSoma(kk,2)/res, 80-(CellSoma(kk,3))/res,'Marker','o', 'MarkerFaceColor','w');
    
    %plot XY
    tempXY = max(volPost,[],3);
    h2 = imagesc(0,90,imrotate(tempXY,-90));set(gca,'YDir','normal'); colormap(jet);
    title(cellIDs{kk});
    plot(160-CellSoma(kk,2)/res, 90 + CellSoma(kk,1)/res,'Marker','o', 'MarkerFaceColor','w');

    %plot yz
    tempYZ = max(volPost,[],1);
    for i = 1:size(tempYZ,3)
        C(i,:) = tempYZ(:,:,i);
    end
    imagesc(160,90,imrotate(C,-90));set(gca,'YDir','normal'); colormap(jet);
    title(cellIDs{kk});
    plot(160+(CellSoma(kk,3))/res, 90+ CellSoma(kk,1)/res,'Marker','o', 'MarkerFaceColor','w');
 
    axis off
    axis image
    
    [m,I] = max(volPost(:));
    maxDesnity(kk) = m;
    [x,y,z] = ind2sub(size(volPost),I);
    XPost(kk) = x; YPost(kk) = y; ZPost(kk) = z-10000/res;
end

hold off;
PostDensityPeak = [YPost'*res,XPost'*res,ZPost'*res];

for i = 1:size(cellIDs,2) 
    lengthToPostPeakNode(i) = findPathLength([cellIDs{i} , '_WithTags.swc'],[5,5,45],PostDensityPeak(i,:),[-1:10]);
end

%%
%Presynapse Heatmap

clear volPre;
clear('tempXY','tempXZ','tempYZ');
clear ('B','C');
figure(15);
for kk = 1:length(allPreSynapse)
    if cellfun('isempty',allPreSynapse(1,kk)) == 1
        continue;
    else
        volPre = zeros((15e4-0)/res,(2.6e5-0)/res,(8e4-0)/res);
        for n = 1:size(allPreSynapse{kk},1);
            volPre(round(allPreSynapse{kk}(n,1)/res) + (-ksize/2:ksize/2), round(allPreSynapse{kk}(n,2)/res) + (-ksize/2:ksize/2), round(allPreSynapse{kk}(n,3)/res + 10000/res)+ (-ksize/2:ksize/2)) = ...
                volPre(round(allPreSynapse{kk}(n,1)/res) + (-ksize/2:ksize/2), round(allPreSynapse{kk}(n,2)/res) + (-ksize/2:ksize/2),round(allPreSynapse{kk}(n,3)/res + 10000/res)+ (-ksize/2:ksize/2))  + K3D;
        end
        subplot(3,8,kk);
        %plot XY
        tempXY = max(volPre,[],3);
        %figure('Units','pixels','Position', [0,0,m,n]);
        imagesc(0,90,imrotate(tempXY,-90));set(gca,'YDir','normal'); colormap(jet);
        title(cellIDs{kk});
        axis off
        hold on;
        plot(160-CellSoma(kk,2)/res, 90 + CellSoma(kk,1)/res,'Marker','o', 'MarkerFaceColor','w');
        axis image;
        %plot XZ
        tempXZ = max(volPre,[],2);
        for i = 1:size(tempXZ,3)
            B(:,i) = tempXZ(:,:,i);
        end
        imagesc(0,0,imrotate(B,-90));set(gca,'YDir','normal'); colormap(jet);
        title(cellIDs{kk});
        axis off
        hold on;
        plot(160-CellSoma(kk,2)/res,80-(CellSoma(kk,3))/res,'Marker','o', 'MarkerFaceColor','w');
        axis image;
        clear B;
        %plot yz
        tempYZ = max(volPre,[],1);
        for i = 1:size(tempYZ,3)
            C(i,:) = tempYZ(:,:,i);
        end
        imagesc(160,90,imrotate(C,-90));set(gca,'YDir','normal'); colormap(jet);
        title(cellIDs{kk});
        axis off
        hold on;
        plot(160+(CellSoma(kk,3))/res, 90+ CellSoma(kk,1)/res,'Marker','o', 'MarkerFaceColor','w');
        axis image;
        clear C;
        
    end
    
    [m,I] = max(volPre(:));
    maxPreDesnity(kk) = m;
    [x,y,z] = ind2sub(size(volPre),I);
    XPre(kk) = x; YPre(kk) = y; ZPre(kk) = z-10000/res;
      
end
PreDensityPeak = [YPre'*res,XPre'*res,ZPre'*res];

for i = 1:size(cellIDs,2)
    if cellfun('isempty',allPreSynapse(1,i)) == 1
        continue;
    else
        lengthToPrePeakNode(i) = findPathLength([cellIDs{i} , '_WithTags.swc'],[5,5,45],PreDensityPeak(i,:),[-1:10]);
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
   % text(lengthToPostPeakNode(i),rho(i)+0.05,cellIDs{i});
end
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



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Distribution of Synapse onto Cells
%figure();
%[m,Index] = sort(rho);
%cellIDsSorted = cellIDs(Index);
%allPostSorted = allPost(Index);
%allRawLengthSorted = allRawLength(Index);
% distribution of postsynaptic sites
figure(18);
for i = 1:length(cellIDs)
    lengthToPostNode = findPathLength([cellIDs{i} , '_WithTags.swc'],[5,5,45],allPost{i},[-1:10]);
    allLengthToPostNode{i} = lengthToPostNode ;
    subplot(3,8,i);
    scatter(1:size(allPost{i},1),sort(lengthToPostNode)/allRawLength(i));
    title(cellIDs{i})
end
 figure(19);
 for i = 1:length(cellIDs)
 subplot(3,8,i);
 hist(sort(allLengthToPostNode{i})/allRawLength(i), length(allPost{i}));title(cellIDs{i});
 
 end
 
 % distribution of presynaptic sites
 figure(20);
 for i = 1:length(cellIDs)
      if cellfun('isempty',allPreSynapse(1,i)) == 1
        continue;
    else
    lengthToPreNode = findPathLength([cellIDs{i} , '_WithTags.swc'],[5,5,45],allPreSynapse{i},[-1:10]);
    allLengthToPreNode{i} = lengthToPreNode ;
    subplot(3,8,i);
    scatter(1:size(allPreSynapse{i},1),sort(lengthToPreNode)/allRawLength(i));
    title(cellIDs{i})
      end
end
 figure(21);
 for i = 1:length(cellIDs)
     if cellfun('isempty',allPreSynapse(1,i)) == 1
        continue;
    else
 subplot(3,8,i);
 hist(sort(allLengthToPreNode{i})/allRawLength(i), length(allPreSynapse{i}));
 title(cellIDs{i});
     end
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Misc. plots
% plot emperical CDFs for all cells
for i = 1:22
p = ((1:size(allLengthToPostNode{i},1))-0.5)'./size(allLengthToPostNode{i},1);
stairs(sort(allLengthToPostNode{i})/allRawLength(i),p);
hold on;
end

% plot rawLength vs. number of PostSynapses

for i = 1:length(cellIDs)
    if ismember(cellIDs{i},cellIDsAlx)==1
        plot(allRawLength(i),length(allPost{i}),'*','Color',[1 0.5 0]);
    elseif ismember(cellIDs{i},cellIDsDbx)==1
        plot(allRawLength(i),length(allPost{i}),'*','Color',[1 0 1]);
    else
        plot(allRawLength(i),length(allPost{i}),'*','Color',[0 0.5 1]);
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
        plot(i,mean(allLengthToPostNode{i}/allRawLength(i)) ,'*','Color',[1 0.5 0]);
        
    elseif ismember(cellIDs{i},cellIDsDbx)==1
        plot(i,mean(allLengthToPostNode{i}/allRawLength(i)) ,'*','Color',[1 0 1]);
    else
        plot(i,mean(allLengthToPostNode{i}/allRawLength(i)) ,'*','Color',[0 0.5 1]);
    end
    hold on;
    %text(allRawLength(i),length(allPost{i}),cellIDs{i});
end



