%% plot spatial and temporal neighbor pairs

map = colormap(parula(length(cellIDs)));

% Pairwise Euclidean distance between cells
PairEucDist= squareform(pdist(CellSoma));
[PairEucDistSort,ICell]  = sort(PairEucDist,2);

%plot all nearest neighbor pairs
for i = 1:22
subplot(3,8,i);
[SomaDistNearCell(i), rhoDistNearCell(i)] = PlotPairs(allTrees{ICell(i,1)}, allTrees{ICell(i,2)}, rho(ICell(i,1)),rho(ICell(i,2)),rho, map);
title(sprintf('%s and %s \n rho1: %.2f and rho2: %.2f', cellIDs{ICell(i,1)}, cellIDs{ICell(i,2)}, rho(ICell(i,1)), rho(ICell(i,2)) ));
end

%axis vis3d;
set(gcf,'color','w');
set(gca,'CLim',[min(rho) max(rho)]);
h = colorbar('eastoutside');
h.Position = [0.7323 0.1121 0.0052 0.2045];
figtitle('Nearest spatial Pair');


% Cells with similar time constants
figure;

PairTimeDist = squareform(pdist(rho'));
[PairTimeDistSort,IRho]  = sort(PairTimeDist,2);

for i = 1:22
subplot(3,8,i);
[SomaDistNearTime(i), rhoDistNearTime(i)] = PlotPairs(allTrees{IRho(i,1)}, allTrees{IRho(i,2)}, rho(IRho(i,1)),rho(IRho(i,2)),rho, map);
title(sprintf('%s and %s \n rho1: %.2f and rho2: %.2f', cellIDs{IRho(i,1)}, cellIDs{IRho(i,2)}, rho(IRho(i,1)), rho(IRho(i,2)) ));
end


%axis vis3d;
set(gcf,'color','w');
set(gca,'CLim',[min(rho) max(rho)]);
h=colorbar('eastoutside');
h.Position = [0.7323 0.1121 0.0052 0.2045];
figtitle('Nearest temporal pair');

%coorelation

figure;
X = [ones(length(SomaDistNearCell),1) SomaDistNearCell'];
y = SomaDistNearTime';
b = X\y;
yCalc2 = X*b;
Rsq2 = 1 - sum((y - yCalc2).^2)/sum((y - mean(y)).^2);

plot(SomaDistNearCell, SomaDistNearTime,'ob', 'MarkerFaceColor' ,'b');
xlabel('Nearest spatial cell (nm)');
ylabel('Nearest temporal cell (nm)');
hold on;
plot(SomaDistNearCell',yCalc2,'-r');
text(max(SomaDistNearCell'),min(yCalc2), 'R^2=0.04');
set(gcf,'color','w');

%% 

map = colormap(parula(length(cellIDs)));

PairEucDist= squareform(pdist(CellSoma));
[PairEucDistSort,ICell]  = sort(PairEucDist,2);

for i = 1:22
subplot(3,8,i);
[SomaDistNearCell(i) tauDistNearCell(i)] = PlotPairs(allTrees{ICell(i,1)}, allTrees{ICell(i,2)}, tau(ICell(i,1)),tau(ICell(i,2)),tau, map);
title(sprintf('%s and %s \n tau1: %.2f and tau2: %.2f', cellIDs{ICell(i,1)}, cellIDs{ICell(i,2)}, tau(ICell(i,1)), tau(ICell(i,2))),'FontName','Arial','FontWeight','normal');
end
h = colorbar;
h.Limits = [min(tau) max(tau)];
h.Location = 'manual';
h.Position = [0.7323 0.1121 0.0052 0.2045];
set(gcf,'color','w');
figtitle('Nearest spatial Pair');

figure;
PairTimeDist = squareform(pdist(tau'));
[PairTimeDistSort,ITau]  = sort(PairTimeDist,2);
for i = 1:22
subplot(3,8,i);
[SomaDistNearTime(i) tauDistNearTime(i)] = PlotPairs(allTrees{ITau(i,1)}, allTrees{ITau(i,2)}, tau(ITau(i,1)),tau(ITau(i,2)),tau, map);
title(sprintf('%s and %s \n tau1: %.2f and tau2: %.2f', cellIDs{ITau(i,1)}, cellIDs{ITau(i,2)}, tau(ITau(i,1)), tau(ITau(i,2))),'FontName','Arial','FontWeight','normal');
end

h = colorbar;
h.Limits = [min(tau) max(tau)];
h.Location = 'manual';
h.Position = [0.7323 0.1121 0.0052 0.2045];
figtitle('Nearest temporal pair');
set(gcf,'color','w');


figure;
%% Plot pairwise Eucledian distances between somas vs Pairwise difference in persistence meausre
clear RhoDiff; clear RhoDiffSD;
clear MeanPairEucDist;
clear RhoRatio;
clear RhoRatioSD

PairEucDist= squareform(pdist(CellSoma));

index = 1;
steps = 1:10000:max(PairEucDist(:));
for i = 1:length(steps)-1
    temp = find(PairEucDist>steps(i) & PairEucDist<steps(i+1));
    [m,n] = ind2sub(size(PairEucDist),temp);
    RhoDiff(index) = mean(abs(rho(m)-rho(n)));
    tempdiff = abs(rho(m)-rho(n));
    RhoRatio(index) = mean(rho(m)./rho(n));
    tempratio = rho(m)./rho(n);
    RhoDiffSD(index) = std(abs(rho(m)-rho(n)));
    RhoRatioSD(index) = std(rho(m)./rho(n));
    MeanPairEucDist(index) = mean(PairEucDist(temp));
    index = index+1;
    figure(1);
    plot(MeanPairEucDist(i),tempdiff,'.r');
    hold on;
    figure(2);
    plot(MeanPairEucDist(i), tempratio, '.r');
    hold on;
end
figure(1);
%plot(MeanPairEucDist,RhoDiff,'o');
errorbar(MeanPairEucDist,RhoDiff,RhoDiffSD,'ok','MarkerFaceColor','k');

xlabel('Pairwise Eucledian distance in nm');
ylabel('Pairwise difference in persistence measure \rho');
clear temp;

figure(2);
errorbar(MeanPairEucDist,RhoRatio,RhoRatioSD,'ok','MarkerFaceColor','k')
xlabel('Pairwise Eucledian distance in nm');
ylabel('Pairwise persistence ratio');

%% Pair wise Eucledian distance in DV and RC axis
% for RC axis
temp = [CellSoma(:,1),CellSoma(:,2)]; % considering only the x,y coordinates
RCpdist = tril(squareform(pdist(temp)),-1);
clear RhoDiff;
clear RhoDiffSD;
clear MeanRCEucDist;

index = 1;
steps = 1:10000:max(RCpdist(:));
for i = 1:length(steps)-1
    temp = find(RCpdist>steps(i) & RCpdist<steps(i+1));
    [m,n] = ind2sub(size(RCpdist),temp);
    %RhoDiff(index,1:length(temp)) = abs(rho(m)-rho(n));
    RhoDiff(index) = mean(abs(rho(m)-rho(n)));
    RhoDiffSD(index) = std(abs(rho(m)-rho(n)));
    %MeanPairEucDist(index,1:length(temp)) = PairEucDist(temp);
    MeanRCEucDist(index) = mean(RCpdist(temp));
    index = index+1;
end
figure;
%plot(MeanPairEucDist,RhoDiff,'o');
errorbar(MeanRCEucDist,RhoDiff,RhoDiffSD,'o')
xlabel('Pairwise RC Eucledian distance in nm');
ylabel('Pairwise difference in persistence measure \rho');
clear temp;

%%
% for DV axis

clear MeanDVEucDist;
clear RhoDiff;
clear RhoDiffSD;
DVpdist = tril(squareform(pdist(CellSoma(:,3))),-1);

index = 1;
steps = 1:2500:max(DVpdist(:));
for i = 1:length(steps)-1
    temp = find(DVpdist>steps(i) & DVpdist<steps(i+1));
    [m,n] = ind2sub(size(DVpdist),temp);
    %RhoDiff(index,1:length(temp)) = abs(rho(m)-rho(n));
    RhoDiff(index) = mean(abs(rho(m)-rho(n)));
    RhoDiffSD(index) = std(abs(rho(m)-rho(n)));
    %MeanPairEucDist(index,1:length(temp)) = PairEucDist(temp);
    MeanDVEucDist(index) = mean(DVpdist(temp));
    index = index+1;
end
figure;
%plot(MeanPairEucDist,RhoDiff,'o');
errorbar(MeanDVEucDist,RhoDiff,RhoDiffSD,'o')
xlabel('Pairwise DV Eucledian distance in nm');
ylabel('Pairwise difference in persistence measure \rho');
clear temp;


%% Plot Pairwise Eucledian distance between somas vs Pairwise difference in number of post synaptic sites
figure();
PairEucDist= tril(squareform(pdist(CellSoma)));
clear temp2;
index = 1;
steps = 1:10000:max(PairEucDist(:));
for i = 1:length(steps)-1
    temp = find(PairEucDist>steps(i) & PairEucDist<steps(i+1));
    [m,n] = ind2sub(size(PairEucDist),temp);
    for j = 1:length(m)
        temp2(j) = abs(length(allPost{m(j)})-length(allPost{n(j)}));
        %temp2(j) = mean(length(allPost{m(j)})+ length(allPost{n(j)}));
    end
    %SynDiffSE(index) = mean(temp2)/sqrt(length(temp2));
    SynDiffSE(index) = std(temp2);
    SynDiffMean(index) = mean(temp2);
    MeanPairEucDist(index) = mean(PairEucDist(temp));
    index = index+1;
    %     plot(MeanPairEucDist(i),temp2,'.r');
    %     hold on;
    clear temp2;
end
%plot(MeanPairEucDist,SynDiff,'o');
errorbar(MeanPairEucDist,SynDiffMean,SynDiffSE,'o');
xlabel('Pairwise Eucledian distance in nm');
ylabel('Pairwise difference in postsynaptic sites ');
%% plot pairwise euclidean distance for each axis (DV, RC)

% RC axis
clear MeanPairRCEucDist;
clear SynDiffMean;
clear SynDiffSE;
figure();

index = 1;
steps = 1:10000:max(RCpdist(:));
for i = 1:length(steps)-1
    temp = find(RCpdist>steps(i) & RCpdist<steps(i+1));
    [m,n] = ind2sub(size(RCpdist),temp);
    for j = 1:length(m)
        temp2(j) = abs(length(allPost{m(j)})-length(allPost{n(j)}));
        %temp2(j) = mean(length(allPost{m(j)})+ length(allPost{n(j)}));
    end
    %SynDiffSE(index) = mean(temp2)/sqrt(length(temp2));
    SynDiffSE(index) = std(temp2);
    SynDiffMean(index) = mean(temp2);
    MeanPairRCEucDist(index) = mean(RCpdist(temp));
    index = index+1;
    %     plot(MeanPairRCEucDist(i),temp2,'.r');
    %     hold on;
    clear temp2;
end
%plot(MeanPairEucDist,SynDiff,'o');
errorbar(MeanPairRCEucDist,SynDiffMean,SynDiffSE,'o');
xlabel('Pairwise RC Eucledian distance in nm');
ylabel('Pairwise difference in postsynaptic sites ');

%% plot pairwise euclidean distance for each axis (DV, RC)

% DV axis
clear MeanPairDVEucDist;
clear SynDiffMean;
clear SynDiffSE;
figure();

index = 1;
steps = 1:2500:max(DVpdist(:));
for i = 1:length(steps)-1
    temp = find(DVpdist>steps(i) & DVpdist<steps(i+1));
    [m,n] = ind2sub(size(DVpdist),temp);
    for j = 1:length(m)
        temp2(j) = abs(length(allPost{m(j)})-length(allPost{n(j)}));
        %temp2(j) = mean(length(allPost{m(j)})+ length(allPost{n(j)}));
    end
    %SynDiffSE(index) = mean(temp2)/sqrt(length(temp2));
    SynDiffSE(index) = std(temp2);
    SynDiffMean(index) = mean(temp2);
    MeanPairDVEucDist(index) = mean(DVpdist(temp));
    index = index+1;
    %     plot(MeanPairDVEucDist(i),temp2,'.r');
    %     hold on;
    clear temp2;
end
%plot(MeanPairEucDist,SynDiff,'o');
errorbar(MeanPairDVEucDist,SynDiffMean,SynDiffSE,'o');
xlabel('Pairwise DV Eucledian distance in nm');
ylabel('Pairwise difference in postsynaptic sites ');

%% plot convexhull of a tree and display, does not display cell

tree = allTrees{21};
for kk = 1:numel(tree)
    children = tree{kk}{2};
    for mm=1:numel(children)
        tempx=[tree{children(mm)}{3}(1); tree{children(mm)}{4}{1}(:,1); tree{kk}{3}(1)];
        tempy=[tree{children(mm)}{3}(2); tree{children(mm)}{4}{1}(:,2); tree{kk}{3}(2)];
        tempz=-[tree{children(mm)}{3}(3); tree{children(mm)}{4}{1}(:,3); tree{kk}{3}(3)];
    end
    plot3(tempx,tempy,tempz);
    X = tempx; Y = tempy; Z = tempz;
end
hold on;
ConHull = convhull(X, Y, Z);
trisurf(ConHull,X,Y,Z, 'faceColor','Cyan','FaceAlpha',0.1);


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

% plot heat map of Postsynaptic sites
figure
for kk = 1:length(cellIDs)
    volPost = zeros((15e4-0)/res,(2.6e5-0)/res,(8e4-0)/res);
    % add 10000/res pixels to the zaxis to avoid edge effects
    for n = 1:length(allPost{kk});
        volPost(round(allPost{kk}(n,1)/res) + (-ksize/2:ksize/2), round(allPost{kk}(n,2)/res) + (-ksize/2:ksize/2), round(allPost{kk}(n,3)/res + 10000/res)+ (-ksize/2:ksize/2)) = ...
            volPost(round(allPost{kk}(n,1)/res) + (-ksize/2:ksize/2), round(allPost{kk}(n,2)/res) + (-ksize/2:ksize/2),round(allPost{kk}(n,3)/res + 10000/res)+ (-ksize/2:ksize/2))  + K3D;
    end
    subplot(3,8,kk);
    threeView(volPost,jet);
    title(cellIDs{kk});
    
    %plot XZ soma location
    plot(160-(CellSoma(kk,2))/res, 80 -(CellSoma(kk,3))/res,'Marker','o', 'MarkerFaceColor','w');
    %plot XY soma location
    plot(160-(CellSoma(kk,2))/res, 90 + (CellSoma(kk,1))/res,'Marker','o', 'MarkerFaceColor','w');
    %plot YZ soma location
    plot(160+(CellSoma(kk,3))/res, 90 + (CellSoma(kk,1))/res,'Marker','o', 'MarkerFaceColor','w');
    
    axis off
    axis vis3d
    
    
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

for i = 1:length(cellIDs)
    lengthToPostNode = findPathLength([cellIDs{i} , '_WithTags.swc'],[5,5,45],allPost{i},[-1:10]);
    allLengthToPostNode{i} = lengthToPostNode ;
    subplot(3,8,i);
    scatter(1:size(allPost{i},1),(lengthToPostNode)/allRawLength{i});
    title(cellIDs{i})
end

%% plot the adjacency Matrix for all trees

for i=1:numel(allTrees)
    adjMat{i} = AdjMat(allTrees{i});
end

for i=1:numel(allTrees)
    subplot(3,8,i)
    imagesc(adjMat{i})
    axis square
end

%% Plots the number of treeBranches for each cell

figure();
for i = 1:length(cellIDs)
    hold on
    if ismember(cellIDs{i},cellIDsAlx)==1
        plot(sum(TreeBranches(allTrees{i})),rho(i),'Marker', '*','Color',[1 0.5 0]);
    elseif ismember(cellIDs{i},cellIDsDbx)==1
        plot(sum(TreeBranches(allTrees{i})),rho(i),'Marker', '*','Color',[1 0 1]);
    else
        plot(sum(TreeBranches(allTrees{i})),rho(i),'Marker', '*','Color',[0 0.5 1]);
    end
end
%line([0,300], [0,1]);
title('Rho vs. Number of Branches');
ylabel('Persistence measure Rho' );
xlabel('Number of Branches');

%% Number of terminal nodes for each cell

for ii = 1:length(cellIDs)
    [B,T] = TreeBranches(allTrees{ii});
    T = find(T);
    TerminalNodes = zeros(length(T),3);
    for i = 1:length(T)
        TerminalNodes(i,:) = allTrees{ii}{1,T(i)}{1,3};
    end
    termPathLength{ii} = findPathLength([cellIDs{ii} , '_WithTags.swc'],[5,5,45],TerminalNodes);
end

%% What % of of dendrite is where

figure();
radius = 10000;

for i = 1:size(cellIDs,2)
    tree{i} = load_tree([cellIDs{i} , '_WithTags.swc']);
    shollIntersections{i}= sholl_tree(tree{i},radius); % concentric circles every radius nm
    EucDistPost{i} = pdist2(allTrees{i}{1}{3},allPost{i}(:,:))';
    
    [B,T] = TreeBranches(allTrees{i});
    B = find(B);
    for kk = 1:length(B)
        BranchPoints(kk,:) = allTrees{i}{B(kk)}{3};
    end
    EucDistBranch{i} = pdist2(allTrees{i}{1}{3},BranchPoints(:,:))';
    
    for jj = 1:length(shollIntersections)
        [Y,I] = find(EucDistPost{i}>(jj-1)*radius & EucDistPost{i}<(jj)*radius);
        count(jj+1) = length(I);
        
        [Y2,I2] = find(EucDistBranch{i}>(jj-1)*radius & EucDistBranch{i}<(jj)*radius);
        count2(jj+1) = length(I2);
    end
%     figure(1)
%     subplot(3,8,i);
%     title(cellIDs{i});
%     plot(1:radius/1000:radius*length(shollIntersections) / 1000,shollIntersections/sum(shollIntersections),'-g*', 1:radius/1000:radius*length(count) /1000,count/sum(count), '-r*');
%     legend('ShollInteresctions', 'PostSynapticSites');
%     
%     figure(2)
%     %title('Dendritic Branching % vs Number of PostSynaptic Sites')
%     subplot(3,8,i);
%     title(cellIDs{i});
%     plot(1:radius/1000:radius*length(count2) /1000,count2/sum(count2), '-g*',1:radius/1000:radius*length(count) /1000,count/sum(count), '-r*' );
%     legend('Dendrite','PostSynapticSites');
    
end


%% find intersecting nodes; 

for i = 1:length(cellIDs)
    lin{i} = lineage(allTrees{i},1,1);
    for ii = 1:length(lin{i})
        IntersectionNodes(ii,:) = allTrees{i}{lin{i}(ii)}{1,3};
    end
    plength{i} = findPathLength([cellIDs{i} , '_WithTags.swc'], [5 5 45], IntersectionNodes);
    plength{i} = plength{i}/1000; % convert from nm to um
end

%% Distrubution of path length of all nodes and distrubution of pathlength of all post synapses

for i = 1:length(cellIDs)
    subplot(3,8,i)
    histogram(plength{i},'BinWidth', 50, 'Normalization','probability');
    hold on;
    histogram(allLengthToPostNode{i}/1000,'BinWidth', 50, 'Normalization','probability');
end




