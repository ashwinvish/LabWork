%figure();

% %B = 12;
% A = [4,5,6,7,16,21,22];
% clear temp1;
% clear temp2;
% for i = 1:size(A,2)
%     temp1 = allPreSynapse{A(i)}; temp2 = allPostSynapse{A(i)};
%     %subplot(3,3,i);
%     treeVisualizer(allTrees{A(i)}, [1],[eval([cellIDs{A(i)},'_axon'])],[{temp2} {temp1}],true,{[1,0.5,0] [1 0.4 0.4] }, 1:numel(allTrees{A(i)}), false);
%     %temp3 = allPreSynapse{B}; temp4 = allPostSynapse{B};
%     % treeVisualizer(allTrees{B}, [1],[eval([cellIDs{B},'_axon'])],[{temp4} {temp3}],false,{[1,0.5,0]}, 1:numel(allTrees{B}), false);
%     text(0.5,0.98,cellIDs{A(i)},'Units','normalized');
%     export_fig(sprintf('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/ZfishFigures/%s_WithTags.eps', cellIDs{A(i)}));
% end
A = [22]%,5,6,7,16,21,22];
temp1 = allPreSynapse{A}; temp2 = allPostSynapse{A};
treeVisualizer(allTrees{A}, [1],[eval([cellIDs{A},'_axon'])],[{temp2} {temp1}],true,{[1,0.5,0] [1 0.4 0.4] }, 1:numel(allTrees{A}), false);
hold on;
text(0.5,0.98,cellIDs{A},'Units','normalized');

for i = 1:length(cellIDs)
    for j = 1:length(allPost{1,i})
        plot3(allPost{1,i}(j,1),allPost{1,i}(j,2),-allPost{1,i}(j,3), 'Marker','o', 'MarkerSize' , 3, 'LineWidth', 0.1 , 'MarkerFaceColor', [0 0.4 0.4] , 'MarkerEdgeColor' , 'none' );
    end
end
view (-140,18);

export_fig(sprintf('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/ZfishFigures/%s_WithAllSynapses.eps', cellIDs{A}));

% for kk = 1:numel(cellIDs)
%     subplot(3,8,kk);
%     if ismember(cellIDs{kk},cellIDsDbx)==1
%         treeVisualizer(allTrees{kk}, [1],[eval([cellIDs{kk},'_axon'])],[allPost(kk) allPreSynapse(kk)],false,{[1, 0, 1] [0.6 0 1]}, 1:numel(allTrees{kk}), false); % Magenta for contraxon
%         treeVisualizer(thisTree, [1],[],[],false,{[1, 0, 1]}, 1:numel(thisTree), false); % Magenta for Dbx
%     else
%         ismember(cellIDs{kk},cellIDsL)==1
%         treeVisualizer(allTrees{kk}, [1],[eval([cellIDs{kk},'_axon'])],[allPost(kk) allPreSynapse(kk)],false,{[ 0.4588 0.7333 0.9922] [0.2941 0.3647 0.0863]}, 1:numel(allTrees{kk}), false); % Blue for neither
%         treeVisualizer(thisTree, [1],[],[],false,{[0, 0.5, 1]}, 1:numel(thisTree), false); % Blue for lateral
%     end
%     text(0.5,0.98,cellIDs{kk},'Units','normalized');
%     
% end
% axis vis3d;
% h1 = gcf;
% h2 = PlotViews(h1);

%export_fig(sprintf('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/ZfishFigures/%s_WithTags.eps', cellIDs{kk}));

%%
%for i=1:numel(allTrees)
%cellIDsAlx = {'Int1_4','Int1_5','Int1_6','Int1_7','Int2_6','Int2_9','Int3_5','Int3_6'};
%cellIDsAlx = [4,5,6,7,13,19,21,22];


% for i = 1:length(cellIDs)
%   if  ismember(cellIDs{i},cellIDsAlx)==1
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
axis square;

for i = 1:size(cellIDs,2)
    overlapNormalized (i,:) = overlap(i,:)/max(overlap(i,:));
end
figure; imagesc(overlapNormalized); axis square;

%%

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
        figtitle('Synaptic pairs <5\mum apart');
        clear a;
        clear b;
    else
        continue
    end
end
export_fig('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/ZfishFigures/PairwiseSynapseDistrubution', '-eps');

% plot the minimum synaptic distance between trees
figure();
imagesc(minDistance/1000); % in microns
colorbar;
title('Minimum synaptic distance between trees in \mum');
xlabel('PostSynaptic cell');
ylabel('PreSynaptic cell');
axis square;
export_fig('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/ZfishFigures/MinimumSynapticDistance','-eps');

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
export_fig('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/ZfishFigures/MinimumSynapticDistance1000', '-eps');




