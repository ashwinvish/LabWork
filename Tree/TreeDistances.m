
% figure properties

%figure('units','normalized','outerposition',[0 0 1 1]);
hold on;

% Presynaptic Cell

 A = [13]; %

 Axntree = allTrees{A};
% 
% % Plot the presynaptic cell
% 
%  temp1 = allPreSynapse{A}; temp2 = allPostSynapse{A};
  %treeVisualizer(allTrees{A}, [1],[eval([cellIDs{A},'_axon'])],[],false,{[1 0.5 0] [1 0 0] }, 1:numel(allTrees{A}), false);

% consider only axon of the tree

inducingNodes = eval([cellIDs{A},'_axon']);
axonTree = [];

for kk=1:numel(Axntree)
    children = Axntree{kk}{2};
    for mm = 1:numel(children)
        if ismember(kk,inducingNodes)
            Axtempx=[Axntree{children(mm)}{3}(1); Axntree{children(mm)}{4}{1}(:,1); Axntree{kk}{3}(1)];
            Axtempy=[Axntree{children(mm)}{3}(2); Axntree{children(mm)}{4}{1}(:,2); Axntree{kk}{3}(2)];
            Axtempz=-[Axntree{children(mm)}{3}(3); Axntree{children(mm)}{4}{1}(:,3); Axntree{kk}{3}(3)];
            h1 = plot3(Axtempx,Axtempy,Axtempz,'color',[1 0 0],'lineWidth',1);
            %h1.Color(4) = 0.2; %  transparency 0-1
            axonTree = [axonTree; Axtempx Axtempy Axtempz];
        else
            continue;
        end
    end
end

% consider dendritic trees of all remaining cells

remTrees = 1:numel(cellIDs);
remTrees(A) = [];
denTree = [];


for ii = 1:numel(remTrees)
    DenTree = allTrees{remTrees(ii)};
    validNodes = [1:numel(DenTree)];
    for jj = 1:numel(DenTree)
        children = DenTree{jj}{2};
        for nn = 1:numel(children)
            DnTempx = [DenTree{children(nn)}{3}(1); DenTree{children(nn)}{4}{1}(:,1)];% DenTree{nn}{3}(1)];
            DnTempy = [DenTree{children(nn)}{3}(2); DenTree{children(nn)}{4}{1}(:,2)];% DenTree{nn}{3}(2)];
            DnTempz = -[DenTree{children(nn)}{3}(3); DenTree{children(nn)}{4}{1}(:,3)];% DenTree{nn}{3}(3)];
            if ismember(jj,validNodes)
                h2 = plot3(DnTempx,DnTempy,DnTempz,'-','color',[0.8 0.8 0.8]);%,'lineWidth',0.5);
                h2.Color(4) = 0.2;
            end
            denTree = [denTree; DnTempx DnTempy DnTempz];
        end
    end
end


box on;
axis([ 20000 140000 60000 250000 -60000 0]);
plot( [20000, 40000], [70000, 70000],'-k' ) % insert 20um sclaebar
daspect([1 1 1]); % make aspect ratio [1 1 1]

% calculate euclidean distance between axTree and denTree

clear Dist;
distThreshold = 1000; % set threshold for distance
Dist = pdist2(axonTree, denTree);
[r,c] = find(Dist>0 & Dist<distThreshold);

% plot only dendrites inside 1000nm
%plot3(denTree(c,1),denTree(c,2),denTree(c,3),'.','color',[0 0.5 1]);

% unique locations on axon

sprintf('Number of locations on the axon:%d \n Number of locations on dendrites: %d', length(r), length(c))
sprintf('Number of unique locations on the axon:%d \n Number of unique locations on dendrites: %d', length(unique(r)), length(unique(c)))
% 
diffUniq = find(diff(unique(c))>1)
uniquec = unique(c);
s1 = scatter3(denTree(uniquec(diffUniq),1),denTree(uniquec(diffUniq),2),denTree(uniquec(diffUniq),3),'Marker','o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor','k');
drawnow;
s1Markers = s1.MarkerHandle; 
s1Markers.FaceColorData = uint8(255*[1;0;0;0.5]);  % Alpha=0.3 => 70% transparent red
s1.SizeData = 20;


% plot actual synapses

syn = [16662	24415	470
18872	23363	468
19645	36619	960
19279	35061	1046];
syn(:,3) = -1*syn(:,3);
s2 = scatter3(5*syn(:,1),5*syn(:,2),45*syn(:,3),'Marker','o','MarkerFaceColor',[0.2 0.2 1],'MarkerEdgeColor','none');
% s2Markers = s2.MarkerHandle; 
% s2Markers.FaceColorData = uint8(255*[1;0;0;0.3]);  % Alpha=0.3 => 70% transparent red
s2.SizeData = 80;




% drawnow;
% hMarkers = hSyn.MarkerHandle; 
% hMarkers.FaceColorData = uint8(255*[1;0;0;0.3]);  % Alpha=0.3 => 70% transparent red

%clear hMarkers;

Dist1000 = [];
for i = 1:length(r)
    Dist1000 = [Dist1000; Dist(r(i),c(i))];
end
%hist(Dist(r,:)/1000);

% set figure to desired orientation

set (gca,'XTick',[], 'YTick',[],'ZTick', [],'Ydir','reverse');
view([-180,90]);
axis vis3d;
%text(0.5,0.98,cellIDs{A},'Units','normalized');





% save figure

%export_fig(sprintf('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/ZfishFigures/%s_With%dnmDendrites.png',cellIDs{A}, distThreshold),'-r300');




