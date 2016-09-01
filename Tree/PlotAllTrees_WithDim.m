%%
clear all;
loadTrees; % Load all data

%% Display Control Cells

DisplayTree(TreeC1,[1],true);
DisplayTree(TreeC2,[1],false);
DisplayTree(TreeC3,[1],false);
DisplayTree(TreeC4,[1],false);

ControlCellSoma(1,:) =  TreeC1{1,1}{1,3};
ControlCellSoma(2,:) =  TreeC2{1,1}{1,3};
ControlCellSoma(3,:) =  TreeC3{1,1}{1,3};
ControlCellSoma(4,:) =  TreeC4{1,1}{1,3};

hold on;
stripe1 = [ControlCellSoma(1,:);ControlCellSoma(2,:);ControlCellSoma(3,:);ControlCellSoma(4,:)];	%stripe1, corresponds with alx transcription factor
line(stripe1(:,1),stripe1(:,2),-stripe1(:,3),'LineWidth',2,'LineStyle','--');

DisplayTree(TreeC5,[1],false);
DisplayTree(TreeC6,[1],false);
DisplayTree(TreeC7,[1],false);

ControlCellSoma(5,:) =  TreeC5{1,1}{1,3};
ControlCellSoma(6,:) =  TreeC6{1,1}{1,3};
ControlCellSoma(7,:) =  TreeC7{1,1}{1,3};

stripe2 = [ControlCellSoma(7,:);ControlCellSoma(5,:);ControlCellSoma(6,:)];                     	%stripe2, corresponds with dbx transcription factor
line(stripe2(:,1),stripe2(:,2),-stripe2(:,3),'LineWidth',2,'LineStyle','--');

axis vis3d
PlotViews(gcf);

%% Plot all Cells in three axes
for kk = 1:numel(cellIDs)
    %subplot(3,8,kk);
    if ismember(cellIDs{kk},cellIDsAlx)==1
        DisplayTree(allTrees{kk},[1],false,eval([cellIDs{kk},'_axon']),calx, allPreSynapse{kk}, allPost{kk}) % plot with axon hilighted
        %DisplayTree(allTrees{kk},[1],false,eval([cellIDs{kk},'_axon']),[1 0.5 0.3]);     % plot without hilighting axon
        %DisplayTree(allTrees{kk},[1],false,eval([cellIDs{kk},'_axon']),[rand rand rand], allPreSynapse{kk}, allPost{kk}); 
    elseif ismember(cellIDs{kk},cellIDsTrans)==1
        DisplayTree(allTrees{kk},[1],false,eval([cellIDs{kk},'_axon']),ctrans, allPreSynapse{kk}, allPost{kk});
    elseif ismember(cellIDs{kk},cellIDsDbx)==1
        DisplayTree(allTrees{kk},[1],false,eval([cellIDs{kk},'_axon']),cdbx, allPreSynapse{kk}, allPost{kk});
        %DisplayTree(allTrees{kk},[1],false,eval([cellIDs{kk},'_axon']),[1 0.3 1]);
        %DisplayTree(allTrees{kk},[1],false,eval([cellIDs{kk},'_axon']),[rand rand rand], allPreSynapse{kk}, allPost{kk}); 
    else
        DisplayTree(allTrees{kk},[1],false,eval([cellIDs{kk},'_axon']),cbarhl, allPreSynapse{kk}, allPost{kk});
       % DisplayTree(allTrees{kk},[1],false,eval([cellIDs{kk},'_axon']),[0.3 0.5 1]);
       % DisplayTree(allTrees{kk},[1],false,eval([cellIDs{kk},'_axon']),[rand rand rand], allPreSynapse{kk}, allPost{kk}); 

    end
end
hold on;
%scatter3(MauthnerCell(1,1),MauthnerCell(1,2),MauthnerCell(1,3), 500,'p','MarkerFaceColor','k', 'MarkerEdgeColor', 'k'); % location of the Mauther Cell center
%  line(stripe1(:,1),stripe1(:,2),-stripe1(:,3),'LineWidth',2,'LineStyle','-', 'color','k');                  	% stripe1
%  line(stripe2(:,1),stripe2(:,2),-stripe2(:,3),'LineWidth',2,'LineStyle','-','color','k' );                      % stripe2

axis vis3d;
%h1 = gcf;
%h2 = PlotViews(h1);                                                                                             % plots the three views of the figure handle h1

%% Plot all alx cells
index = 0;
clear ha;
%pause(2);
%ha = tight_subplot(2,5,[.05 .05],[.05 .1],[.01 .01]);
for kk = 1:numel(cellIDs)
    if ismember(cellIDs{kk},cellIDsAlx)==1
        index = index+1;
        %axes(ha(index));
        %DisplayTree(allTrees{kk},[1],false,eval([cellIDs{kk},'_axon']),[1 0.5 0.3],allPreSynapse{kk}, allPost{kk});
        DisplayTree(allTrees{kk},[1],false,eval([cellIDs{kk},'_axon']),calx);
       % view(-150,35);
       %x title(sprintf('Cell ID %s', cellIDs{kk}),'FontName','Arial','FontWeight','normal' );
    else
        continue;
    end
end

% set(ha(1:10),'BoxStyle','full');
% set(ha(1:10),'XColor',[0.831, 0.816, 0.784]);
% set(ha(1:10),'YColor',[0.831, 0.816, 0.784]);
% set(ha(1:10),'ZColor',[0.831, 0.816, 0.784]);
% title(sprintf('CellID %s', cellIDs{kk}));
set(gcf,'color','w');
%figtitle('All ipsilateral projecting cells','FontName','Arial', 'FontWeight','Bold');

%% Plot all dbx1b cells
index = 0;
clear ha;
ha = tight_subplot(2,4,[.05 .05],[.05 .1],[.01 .01]);
for kk = 1:numel(cellIDs)
    if ismember(cellIDs{kk},cellIDsDbx)==1
        index = index+1;
        axes(ha(index));
        DisplayTree(allTrees{kk},[],false,eval([cellIDs{kk},'_axon']),cdbx, allPreSynapse{kk}, allPost{kk});
        hold on
        TreeSomata(kk, cdbx);
        %DisplayTree(allTrees{kk},[1].true,[1 0.3 0.1]);
        view(-150,35);
        title(sprintf('Cell ID %s', cellIDs{kk}),'FontName','Arial','FontWeight','normal' );
    else
        continue;
    end
end
set(ha(1:8),'BoxStyle','full');
set(ha(1:8),'XColor',[0.831, 0.816, 0.784]);
set(ha(1:8),'YColor',[0.831, 0.816, 0.784]);
set(ha(1:8),'ZColor',[0.831, 0.816, 0.784]);
title(sprintf('CellID %s', cellIDs{kk}),'FontName','Arial','FontWeight','normal' );
set(gcf,'color','w');
%figtitle('All contralateral projecting cells','FontName','Arial', 'FontWeight','Bold');

%% Plot all barhl1 cells

index = 0;
clear ha;
%ha = tight_subplot(2,3,[.05 .05],[.05 .1],[.01 .01]);
for kk = 1:numel(cellIDs)
    if ismember(cellIDs{kk},cellIDsL)==1
        index = index+1;
 %       axes(ha(index));
        DisplayTree(allTrees{kk},[1],false,eval([cellIDs{kk},'_axon']),cbarhl);
        %DisplayTree(allTrees{kk},[1],true,[0.3 0.5 1]);
 %       view(-150,35);
 %       title(sprintf('Cell ID %s', cellIDs{kk}),'FontName','Arial','FontWeight','normal' );
    else
        continue;
    end
end
% set(ha(1:6),'BoxStyle','full');
% set(ha(1:6),'XColor',[0.831, 0.816, 0.784]);
% set(ha(1:6),'YColor',[0.831, 0.816, 0.784]);
% set(ha(1:6),'ZColor',[0.831, 0.816, 0.784]);
set(gcf,'color','w');
%figtitle('All unknown projecting cells','FontName','Arial', 'FontWeight','Bold');
%% Plot All Trans Cells
index = 0;
clear ha;
%ha = tight_subplot(2,3,[.05 .05],[.05 .1],[.01 .01]);
for kk = 1:numel(cellIDs)
    if ismember(cellIDs{kk},cellIDsTrans)==1
        index = index+1;
 %       axes(ha(index));
        DisplayTree(allTrees{kk},[1],false,eval([cellIDs{kk},'_axon']),ctrans);
        %DisplayTree(allTrees{kk},[1],true,[0.3 0.5 1]);
 %       view(-150,35);
 %       title(sprintf('Cell ID %s', cellIDs{kk}),'FontName','Arial','FontWeight','normal' );
    else
        continue;
    end
end
% set(ha(1:6),'BoxStyle','full');
% set(ha(1:6),'XColor',[0.831, 0.816, 0.784]);
% set(ha(1:6),'YColor',[0.831, 0.816, 0.784]);
% set(ha(1:6),'ZColor',[0.831, 0.816, 0.784]);
set(gcf,'color','w');
%figtitle('All unknown projecting cells','FontName','Arial', 'FontWeight','Bold');

%% Plot all Soma with Time Constants

CellSoma = zeros(length(cellIDs),3);
for kk = 1 : length(allTrees)
    CellSoma(kk,1:3) =  allTrees{1,kk}{1,1}{1,4}{1};
end

%figure();
scatter3(CellSoma(:,1),CellSoma(:,2),-CellSoma(:,3), 500,log2(tau), 'fill', 'Marker','o', 'MarkerEdgeColor', 'k');
%scatter3(CellSoma(:,1),CellSoma(:,2),-CellSoma(:,3), 150, 'red', 'fill', 'Marker','o' ,'LineWidth', 2,'MarkerEdgeColor', 'k');
hold on;
scatter3(MauthnerCell(1,1),MauthnerCell(1,2),MauthnerCell(1,3), 1000,'p','MarkerFaceColor','k', 'MarkerEdgeColor', 'k');
scatter3(Mid2C1(1,1),Mid2C1(1,2),Mid2C1(1,3), 500,'p','MarkerFaceColor','r', 'MarkerEdgeColor', 'k');
scatter3(Mid2C2(1,1),Mid2C2(1,2),Mid2C2(1,3), 500,'p','MarkerFaceColor','r', 'MarkerEdgeColor', 'k');
scatter3(Mid3C1(1,1),Mid3C1(1,2),Mid3C1(1,3), 500,'p','MarkerFaceColor','b', 'MarkerEdgeColor', 'k');
scatter3(Mid3C2(1,1),Mid3C2(1,2),Mid3C2(1,3), 500,'p','MarkerFaceColor','b', 'MarkerEdgeColor', 'k');
scatter3(CAD(1,1), CAD(1,2), CAD(1,3), 500,'p','MarkerFaceColor','g', 'MarkerEdgeColor', 'k');
scatter3(CAV(1,1), CAV(1,2), CAV(1,3), 500,'p','MarkerFaceColor','g', 'MarkerEdgeColor', 'k');
caxis([min(log2(tau)) max(log2(tau))]);
colormap parula;
box on;
axis([ 20000 140000 60000 250000 -60000 0]);
plot( [120000, 140000], [70000, 70000],'-k' ) % insert 20um sclaebar
daspect([1 1 1]); % make aspect ratio [1 1 1]
%set (gca,'XTick',[], 'YTick',[],'ZTick', [], 'Ydir','reverse');
set(gca,'Ydir','reverse');
view([-128,32]); % xy view
set(gca,'BoxStyle','full');
%CT = cbrewer('seq','OrRd',size(CellSoma,1));
%colormap(CT);
%  line(stripe1(:,1),stripe1(:,2),-stripe1(:,3),'LineWidth',2,'LineStyle','--');                       % stripe1
%  line(stripe2(:,1),stripe2(:,2),-stripe2(:,3),'LineWidth',2,'LineStyle','--');                       % stripe2
axis vis3d;

%  figHandle = gcf;
%  PlotViews(figHandle);
%%

figure();
map = colormap(parula(22));

for i = 1:length(cellIDs);
    A = i;
    temp1 = allPreSynapse{A}; temp2 = allPostSynapse{A};
    cellNo = find(rho(i) == sort(rho));
    DisplayTree(allTrees{i},false,[], map(cellNo,:));
end
h = colorbar;
h.Limits = [min(rho) max(rho)];
h.Location = 'manual';
h.Position = [0.6730    0.1100    0.0117    0.8150];
set(gcf,'color','w');

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
plot(cellClass,cellfun(@length,allPost),'o', 'MarkerFaceColor','b');
hold on;
meanCellPost = [mean(IpsiPost(find(IpsiPost))), mean(ContraPost(find(ContraPost))), mean(UnknownPost(find(UnknownPost))) ];
SdCellPost = [std(IpsiPost(find(IpsiPost))), std(ContraPost(find(ContraPost))), std(UnknownPost(find(UnknownPost)))];
errorbar(meanCellPost,SdCellPost,'ok', 'MarkerFaceColor','k');
xlabel('CellGroups');
ylabel('Number of postsynaptic sites');

figure();
plot(cellClass,cell2mat(allRawLength)/1000,'o', 'MarkerFaceColor', 'b');
hold on;
meanCellLength = [mean(IpsiLength(find(IpsiLength)))/1000, mean(ContraLength(find(ContraLength)))/1000, mean(UnknownLength(find(UnknownLength)))/1000 ]; % in um
SdCellLength = [std(IpsiLength(find(IpsiLength))/1000), std(ContraLength(find(ContraLength))/1000), std(UnknownLength(find(UnknownLength))/1000)];
errorbar(meanCellLength,SdCellLength,'ok', 'MarkerFaceColor','k');
xlabel('CellGroups');
ylabel('Neuronal length in \mum');

%% Get cell Somas and plot pairwise distances
% All soma locations

EucDist = squareform(pdist(CellSoma));                                                              % pairwise Euclidean distance between somas

figure();
imagesc(EucDist);                                                                                   % plot pairwise euclidean distance
colormap gray;
axis square;

Links = linkage(pdist(CellSoma));
figure();
dendrogram(Links,'Labels',cellIDs)                                                                  % dendrogram of cell distances based on pairwise Euclidean distance

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

%% Plot heat map of the density of synapses

% Generate 3D gaussian Kernel
clear volPost;
res = 1000;                                                                                          % downsampling factor
ksize = 12000;                                                                                        % size of Kernel in nm

% plot heat map of Postsynaptic sites
figure(14);
for kk = 1:length(cellIDs)
    subplot(3,8,kk);
    [volPost{kk},m, I] = HeatMapFish(ksize, res, allPost{kk},CellSoma(kk,:), cellIDs{kk},true);
    maxDesnity(kk) = m;
    [x,y,z] = ind2sub(size(volPost{kk}),I);
    XPost(kk) = x; YPost(kk) = y; ZPost(kk) = z;                  % check this line again??
end

hold off;
PostDensityPeak = [XPost'*res,YPost'*res,ZPost'*res];

for i = 1:size(cellIDs,2)
    lengthToPostPeakNode(i) = findPathLength([cellIDs{i} , '_WithTags.swc'],allTrees{i},[5,5,45],PostDensityPeak(i,:));
end

%% Presynapse Heatmap

figure(15);
clear volPre
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
    XPre(kk) = x; YPre(kk) = y; ZPre(kk) = z;
    
end
%PreDensityPeak = [YPre'*res,XPre'*res,ZPre'*res];
PreDensityPeak = [XPre'*res,YPre'*res,ZPre'*res];

for i = 1:size(cellIDs,2)
    if cellfun('isempty',allPreSynapse(1,i)) == 1
        continue;
    else
        lengthToPrePeakNode(i) = findPathLength([cellIDs{i} , '_WithTags.swc'],allTrees{i},[5,5,45],PreDensityPeak(i,:));
    end
end

%% calculate the convexhull volume of two intersecting point clouds
overlap = zeros(size(cellIDs,2));
overlapNormalized = zeros(size(cellIDs,2));
for i = 1:size(cellIDs,2)
    
    if ismember(cellIDs{i},cellIDsAlx)==1
        figure();
        title(cellIDs{i});
        
        for j = 1:size(cellIDs,2)
            subplot(3,8,j);
            B = j;
            clear temp3;
            clear temp4;
            clear temp1;
            clear temp2;
            A = i;
            temp1 = allPreSynapse{A}; temp2 = allPostSynapse{A};
            [C1,CV1,x1,y1,z1,htri1] = TreeConvexHull(allTrees{A},[1],[],[{temp2} {temp1}],false,{[rand,rand,rand]},'green',1:numel(allTrees{A}));
            
            temp3 = allPreSynapse{B}; temp4 = allPostSynapse{B};
            [C2,CV2,x2,y2,z2] = TreeConvexHull(allTrees{B},[1],[],[{temp4} {temp3}],false,{[rand,rand,rand]},'red',1:numel(allTrees{B}));
            
            % to calculate if point cloud of B is inside A
            in = inpolyhedron(htri1.Faces,htri1.Vertices,[x2 y2 z2]);
            if size([x2(in) y2(in) z2(in)],1) == 0
                Vint = 0;
                overlap(i,j) = Vint;
            else if ~size([x2(in) y2(in) z2(in)],1) == 0
                    [Cint,Vint] = convhull(x2(in),y2(in),z2(in));
                    hold on;
                    trimesh(Cint,x2(in),y2(in),z2(in),'FaceColor','b','FaceAlpha',0.5,'EdgeColor','black','EdgeAlpha',0.1);
                    overlap(i,j) = Vint/1e+9; % overlap volume in um3
                end
            end
        end
    else
        continue;
    end
end

figure;
imagesc(overlap);
title('Convexhull volume');
axis square;

for i = 1:size(cellIDs,2)
    overlapNormalized (i,:) = overlap(i,:)/max(overlap(i,:));
end
figure; imagesc(overlapNormalized); axis square;
title('Normalized convexhull volume');


%% DotProduct of two intersecting volumes

for i = 1:size(cellIDs,2)
    if ~cellfun('isempty',allPreSynapse(1,i)) == 1					% do only for cells that have presynaptic sites
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

%%

%alx-alx overlapped area
AlxArea = [];
TransArea = [];
AlxDbxArea = [];
AlxBarhlArea = [];
AlxTransArea = [];

for i = 1:numel(cellIDs)
    for j = 1:numel(cellIDs)
        if ismember(cellIDs{i}, cellIDsAlx) ==1 && ismember(cellIDs{j}, cellIDsAlx) ==1
            AlxArea= [AlxArea, IntArea(i,j)];
        else ismember(cellIDs{i}, cellIDsTrans) ==1 && ismember(cellIDs{j}, cellIDsTrans) ==1
            TransArea = [TransArea, IntArea(i,j)];
        end
    end
end

for i = 1:numel(cellIDs)
    for j = 1:numel(cellIDs)
        if ismember(cellIDs{i}, cellIDsAlx) ==1 && ismember(cellIDs{j}, cellIDsDbx) ==1
            AlxDbxArea = [AlxDbxArea, IntArea(i,j)];
        elseif ismember(cellIDs{i}, cellIDsAlx) ==1 && ismember(cellIDs{j}, cellIDsTrans) ==1
            AlxTransArea - [AlxTransArea, IntArea(i,j)];
        else ismember(cellIDs{i}, cellIDsAlx) ==1 && ismember(cellIDs{j}, cellIDsL) ==1
            AlxBarhlArea = [AlxBarhlArea, IntArea(i,j)];
        end
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
    [lengthToPostNode,PostDiff] = findPathLength_new([cellIDs{i} , '_WithTags.swc'],allTrees{i},[5,5,45],allPost{i});
    allLengthToPostNode{i} = lengthToPostNode ;
    allPostDiff{i} = PostDiff;
    subplot(3,8,i);
    scatter(1:size(allPost{i},1),sort(lengthToPostNode)/cell2mat(allRawLength(i)));
    title(cellIDs{i})
end

% distribution of presynaptic sites
figure(20);
for i = 1:length(cellIDs)
    if ~cellfun('isempty',allPreSynapse(1,i)) == 1
        i
        [lengthToPreNode, PreDiff] = findPathLength_new([cellIDs{i} , '_WithTags.swc'],allTrees{i},[5,5,45],allPreSynapse{i});
        allLengthToPreNode{i} = lengthToPreNode ;
        allPreDiff{i} = PreDiff;
        subplot(3,8,i);
        scatter(1:size(allPreSynapse{i},1),sort(lengthToPreNode)/cell2mat(allRawLength(i)));
        %hold all;
        title(cellIDs{i})
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Misc. plots
% plot emperical CDFs for all cells
figure()
for i = 1:22
    p = ((1:size(allLengthToPostNode{i},1))-0.5)'./size(allLengthToPostNode{i},1);
    stairs(sort(allLengthToPostNode{i})/cell2mat(allRawLength(i)),p);
    hold on;
end

title('Emperical CDF for all cells');

% plot rawLength vs. number of PostSynapses
figure();
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
xlabel('Rawlength (nm)');
ylabel('Number of postsynaptic sites');


% plot cellDepth vs. Rho
figure();
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

xlabel('cell depth (nm)');
ylabel('Persistance measure \rho');

%
figure();
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

xlabel(' cell #');
ylabel('ratio of mean(allLengthToPostNode{i}/cell2mat(allRawLength(i))) ');


figure();
[y,I] = sort(CellSoma(:,3));
scatter(1:22,y,50,rho(I),'filled', 'MarkerEdgeColor','k');
set(gca,'YDir','reverse');
xlabel('Cell#');
ylabel('Depth in nm');
title('Persistence measure \rho Vs. cell depth');

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
histogram(allPostSynapticLength/1000,'BinWidth',10); % dimensions in microns
title('Distribution of Postsynaptic pathlength');
xlabel('Postsynaptic pathlenght in \mum');
ylabel('count');
subplot(122);
histogram(allPreSynapticLength/1000,'BinWidth',10); % dimensions in microns
title('Distribution of Presynaptic pathlength');
xlabel('Presynaptic pathlenght in \mum');
ylabel('count');

figure();
h = histogram(allPostSynapticLength/1000,'FaceColor',[0.9,0,0], 'BinWidth', 10); % dimensions in microns
childHandle = get(h,'Children');
set(childHandle,'FaceAlpha',0.7); % 0 = transparent, 1 = opaque.
hold on;
histogram(allPreSynapticLength/1000,'FaceColor',[0,0.8,0], 'BinWidth',10); % dimensions in microns
axis square
box off;
ylabel('count');
xlabel('Pathlength in \mum');
set(gca,'XLim', [0,250],'FontName','Arial', 'FontSize', 40, 'LineWidth',2);

figure();
h1 = histogram(allPostSynapticLength/1000,'FaceColor',[0.9,0,0]); % dimensions in microns
h1.Normalization = 'probability';
h1.BinWidth = 10;
hold on;
h2 = histogram(allPreSynapticLength/1000,'FaceColor',[0,0.8,0]); % dimensions in microns
h2.Normalization = 'probability';
h2.BinWidth = 10;
axis square
box off;
ylabel('Probability');
xlabel('Pathlength in \mum');
set(gca,'XLim', [0,250],'FontName','Arial', 'FontSize', 40, 'LineWidth',2);

%% Distribution of inter-synaptic distance

clear InterPostSynapticDistance;
InterPostSynapticDistance = [];
figure;
%subplot(1,2,1);

for ii = 1:numel(cellIDs)
    subplot(3,8,ii)
    InterPostSynapticDistance = [InterPostSynapticDistance;nonzeros(allPostDiff{ii})];
    histogram(nonzeros(allPostDiff{ii})/1000,'BinWidth',0.5);
end
% histogram(InterPostSynapticDistance/1000, 'BinWidth',0.5);
% title('Inter-postsynaptic distance');
% xlabel('Distance in \mum');
% ylabel('Count');

clear InterPreSynapticDistance;
InterPreSynapticDistance = [];
tempPre =[];
figure();
%subplot(1,2,2);
index =1;

for ii = 1:size(cellIDs,2)
    if ~isempty(allPreDiff{ii})
        ii;
        InterPreSynapticDistance = [InterPreSynapticDistance; nonzeros(allPreDiff{ii})];
    end
end
             
histogram(InterPreSynapticDistance/1000, 'BinWidth', 0.5)
title('Inter-presynaptic distance');
xlabel('Distance in \mum');
ylabel('Count');



%% ratio of dendritic length/ axonal length
axLength = [];
clear temp;

for i = 1:size(cellIDs,2)
    if eval([cellIDs{i},'_axon'])>0
        AxnNodes = eval([cellIDs{i},'_axon']);
                AxnNodes = sort(AxnNodes);
        tempLength = 0;
        for jj = 1:numel(AxnNodes)
            tempLength = tempLength + sum(allTrees{i}{AxnNodes(jj)}{1,4}{1,2});
        end
        axLength = [axLength,tempLength];
        
    else
        axLength = [axLength,0];
        continue;
    end
end
denLength = cell2mat(allRawLength)-axLength;
sprintf('dendrite length / axon length = %d',sum(denLength)/sum(axLength));


%% plot number of synapses per cell

[y,I] = sort(cellfun(@length,allPost));


for i = 1:length(cellIDs)
    if ismember(cellIDs(I(i)),cellIDsAlx) == 1
        BarCMap(i,:) = calx;
    elseif ismember(cellIDs(I(i)),cellIDsTrans) == 1
        BarCMap(i,:) = ctrans;
    elseif ismember(cellIDs(I(i)),cellIDsDbx) == 1
        BarCMap(i,:) = cdbx;
    else
        BarCMap(i,:)= cbarhl;
    end
    
    h = bar(i,y(i));
    set(h, 'FaceColor', BarCMap(i,:));
    hold on;
    
end

xlabel('Neuron #');
ylabel('Number of post synapses');
title('Number of synpapse for population');

figure();
[ax,h1,h2] = plotyy(1:22 ,y, 1:22,(y./denLength(I)) * 1000, 'bar', 'plot');
h1.FaceColor =  [0.7,0.7,0.7];
h2.Color =  'k';
h2.Marker = 'o';
h2.MarkerFaceColor = 'k'
set(ax(1),'xcolor','k');
set(ax(1),'ycolor',[0.7,0.7,0.7]);
set(ax(2),'ycolor','k');
set(gca,'FontName', 'Arial', 'FontSize', 40);
set(ax(1), 'XLim',[0,23],  'LineWidth', 2);
set(ax(2), 'XLim', [0,23],  'LineWidth', 2);
xlabel('Neuron #', 'FontName', 'Arial', 'FontSize', 40);
ylabel(ax(1),'Number of postsynaptic sites', 'FontName', 'Arial', 'FontSize', 40);
ylabel(ax(2),'Synapse density (#/\mum)', 'FontName', 'Arial', 'FontSize', 40);
set(ax(2),'FontName', 'Arial', 'FontSize', 40,  'LineWidth', 2);
box off;
set(gcf,'color','w');
axis (ax(1),'square');
axis(ax(2),'square');

figure();
[ax,h1,h2] = plotyy(1:22 ,y./sort(cellfun(@sum,Branches)-rs), 1:22,(y./denLength(I)) * 1000, 'bar', 'plot');
h1.FaceColor = [0.7,0.7,0.7];
h2.Color =  'k';
h2.Marker = 'o';
h2.MarkerFaceColor = 'k';
set(ax(1),'xcolor','k');
set(ax(1),'ycolor',[0.7,0.7,0.7]);
set(ax(2),'ycolor','k');
xlabel('Neuron #');
ylabel(ax(1),' No. of postsynaptic sites/ No. of branches', 'FontName', 'Arial', 'FontSize', 40);
ylabel(ax(2),'Synapse density (#/\mum)','FontName', 'Arial', 'FontSize', 40);
set(gca,'FontName', 'Arial', 'FontSize', 40);
set(ax(1), 'XLim',[0,23],  'LineWidth', 2);
set(ax(2), 'XLim', [0,23],  'LineWidth', 2);
set(ax(2), 'FontName', 'Arial', 'FontSize', 40,  'LineWidth', 2);
box off;
set(gcf,'color','w');
axis (ax(1),'square');
axis(ax(2),'square');

[y,I] = sort(cellfun(@length,allPreSynapse));
figure()
[ax,h1,h2] = plotyy(1:22 ,y, 1:22,(y./axLength(I)) * 1000,'bar', 'plot');
h1.FaceColor =  [0.7,0.7,0.7];
h2.Color =  'k';
h2.Marker = 'o';
h2.MarkerFaceColor = 'k';
set(ax(1),'xcolor','k');
set(ax(1),'ycolor',[0.7,0.7,0.7]);
set(ax(2),'ycolor','k');
xlabel('Neuron #');
ylabel(ax(1),'Number of presynaptic sites','FontName', 'Arial', 'FontSize', 40);
ylabel(ax(2),'Synapse density (#/\mum)','FontName', 'Arial', 'FontSize', 40);
set(ax(1), 'XLim',[0,23],  'LineWidth', 2);
set(ax(2), 'XLim', [0,23],  'LineWidth', 2);
set(ax(2), 'FontName', 'Arial', 'FontSize', 40,  'LineWidth', 2);
set(gca,'FontName', 'Arial', 'FontSize', 40);
box off;
set(gcf,'color','w');
axis (ax(1),'square');
axis(ax(2),'square');

% figure();
% h = tight_subplot(3,8,[.05 .05],[.05 .1],[.01 .01]);
% for i =1:length(cellIDs)
%     [Branches{i}, Terminals{i}, BranchOrder{i}] = TreeBranches(allTrees{i});
%     axes(h(i));
%     title(sprintf('Cell ID %s', cellIDs{i}));
%     BranchOrderVisualizer(allTrees{i},[1],[BranchOrder{i}]);
% end
% figtitle('Branch Order visualization');
%
% figure();
% h = tight_subplot(3,8,[.05 .05],[.05 .1],[.01 .01]);
% for i =1:length(cellIDs)
%     axes(h(i));
%     histogram(BranchOrder{i},'BinLimits',[min(BranchOrder{i}), max(BranchOrder{i})]);
%     title(sprintf('Cell ID %s', cellIDs{i}));
% end

% figtitle('Branch order distribution for all cells');
%
% figure();
% [rs,cs] = cellfun(@size,allSpine);
% bar(sort(cellfun(@sum,Branches)-rs));
% xlabel('neuron#');
% ylabel('Number of branches per cell');







