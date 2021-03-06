% plot all trees colored with their rho value respectively
map = colormap(parula(22));
%[sortedRho,I] = sort(rho);
index = [];

for i = 1:22;
A = i;
temp1 = allPreSynapse{A}; temp2 = allPostSynapse{A};
index = find(rho(i)==sort(rho));
%treeVisualizer(allTrees{A}, [1],[],[],false,{map(i,:)}, 1:numel(allTrees{A}), false);
%DisplayTree(AxnTree,false,eval([cellIDs{A},'_axon']),map(I(i),:), allPreSynapse{A}, allPostSynapse{A}); % with synapses
DisplayTree(allTrees{A},false,[],map(index,:));                                                          % without synapses
hold on;
%text(0.5,0.98,cellIDs{A},'Units','normalized');
%text(0.1,0.1,cellIDs{A});
end
axis vis3d;
set(gcf,'color','w');
set(gca,'CLim',[min(rho) max(rho)]);
colorbar('eastoutside');
%PlotViews(gcf);
% h = colorbar;
% h.Limits = [min(rho) max(rho)];
% h.Location = 'manual';
% h.Position = [0.6730    0.1100    0.0117    0.8150];


%% plot all ipsi cells and overlaping convex hulls

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
title('Overalap of convexhull volumes');
axis square;

for i = 1:size(cellIDs,2)
    overlapNormalized (i,:) = overlap(i,:)/max(overlap(i,:));
end
figure;
imagesc(overlapNormalized);
axis square;
title('Overlap of normalized convex hull volumes');

%% plot 3 views of convexhull overlaps of any cell pairs

A = 16;
B = 15;
clear temp3;
clear temp4;
clear temp1;
clear temp2;

temp1 = allPreSynapse{A}; temp2 = allPostSynapse{A};
[C1,CV1,x1,y1,z1,htri1] = TreeConvexHull(allTrees{A},[1],[],[{temp2} {temp1}],false,{[rand,rand,rand]},'green',1:numel(allTrees{A}));

temp3 = allPreSynapse{B}; temp4 = allPostSynapse{B};
[C2,CV2,x2,y2,z2] = TreeConvexHull(allTrees{B},[1],[],[{temp4} {temp3}],false,{[rand,rand,rand]},'red',1:numel(allTrees{B}));
in = inpolyhedron(htri1.Faces,htri1.Vertices,[x2 y2 z2]);
if size([x2(in) y2(in) z2(in)],1) == 0
    disp ('no overlap');
else
    [Cint,Vint] = convhull(x2(in),y2(in),z2(in));
    trimesh(Cint,x2(in),y2(in),z2(in),'FaceColor','b','FaceAlpha',0.5,'EdgeColor','black','EdgeAlpha',0.1);
end

h1 = gcf;
PlotViews(h1);


%% dotProduct of two intersecting volumes

for i = 1:size(cellIDs,2)
    if ismember(cellIDs{i},cellIDsAlx)==1
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

%% Minimum distance between synapses between trees

index =0;
k=0;
figure();
for i = 1:size(cellIDs,2)
    if size(allPreSynapse{i},1)>0
        for j = 1:size(cellIDs,2)
            sz = 0;
            [MinDistance{i,j}] = MinSynapticDistance(allPreSynapse{i},allPost{j});
            sz = sz + sum(size(MinDistance{i,j},1)*size(MinDistance{i,j},2));
            minDistance(i,j) = min(min(MinDistance{i,j}));
            PairwiseSynapticDistance(index+1:index+sz,i) = [reshape(MinDistance{i,j},[],1)];
            index = sz;
        end
        clear temp;
        k=k+1;
        subplot(3,3,k);
        [a,b] = find(PairwiseSynapticDistance(:,i)>0 & PairwiseSynapticDistance(:,i)<5000);
        hist(PairwiseSynapticDistance(a,i)./1000); % report in microns
        xlim([0 5]);
        title(cellIDs{i});
        figtitle('Synaptic pairs <5 \mum apart');
        clear a;
        clear b;
    else
        continue
    end
end
%export_fig('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/ZfishFigures/PairwiseSynapseDistrubution', '-eps');

% plot the minimum synaptic distance between trees
figure();
imagesc(minDistance/1000); % in microns
colorbar;
title('Minimum synaptic distance between trees in \mum');
xlabel('PostSynaptic cell');
ylabel('PreSynaptic cell');
axis square;
%export_fig('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/ZfishFigures/MinimumSynapticDistance','-eps');

% plot only those synapses pairs that are within 1000nm of each other
figure();
[r,c] = find(minDistance>0 & minDistance<1000); % find synapses that are atleast 1000nm close to each other
minDistance1000 = zeros(size(cellIDs,2));
for i = 1:size(r,1)
    minDistance1000(r(i),c(i)) = minDistance(r(i),c(i));
end
imagesc(minDistance1000/1000); % in microns
colorbar;
title('Synapses 1\mum apart');
xlabel('PostSynaptic cell');
ylabel('PreSynaptic cell');
axis square;
%export_fig('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/ZfishFigures/MinimumSynapticDistance1000', '-eps');
