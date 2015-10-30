% shuffle axon tree

AxonTree = 5,
Shuffles = 1000;

[A,B] =  UniqueSites(allTrees,cellIDs,AxonTree,false,true);
axonTree = zeros(size(A,1),3,Shuffles);
axonTreeJitter = zeros(size(A,1),3,Shuffles);

for i = 1:1:1000  
[axonTree(1:size(A,1),1:3,i), axonTreeJitter(1:size(A,1),1:3,i)]= UniqueSites(allTrees,cellIDs,AxonTree,false,true);
end


% consider dendritic trees of all remaining cells

remTrees = 1:numel(cellIDs);
remTrees(AxonTree) = [];
denTree = [];
hold on;

for ii = 1:numel(remTrees)
    DenTree = allTrees{remTrees(ii)};                % iterating through all trees
    validNodes =  eval([cellIDs{ii},'_axon']);       % keeping track of axonal nodes
    if isempty(validNodes)                           % if no aoxon check
        validNodes = 1:numel(DenTree);
    end                                             % if no axon then iterate through all nodes
    for jj = 1:numel(DenTree)
        children = DenTree{jj}{2};              % Consider childern of jj node
        for nn = 1:numel(children)              % iterate over all children
            DnTempx = [DenTree{children(nn)}{3}(1); DenTree{children(nn)}{4}{1}(:,1)];% DenTree{nn}{3}(1)];
            DnTempy = [DenTree{children(nn)}{3}(2); DenTree{children(nn)}{4}{1}(:,2)];% DenTree{nn}{3}(2)];
            DnTempz = -[DenTree{children(nn)}{3}(3); DenTree{children(nn)}{4}{1}(:,3)];% DenTree{nn}{3}(3)];
            if length(validNodes)< length(DenTree) && ismember(jj,validNodes)
                continue;
            else
                h2 = plot3(DnTempx,DnTempy,DnTempz,'-','color',[0.8 0.8 0.8]);%,'lineWidth',0.5);
                %h2.Color(4) = 0.2;
            end
            denTree = [denTree; DnTempx DnTempy DnTempz];
        end
    end
end


% % plot the real axon
% plot3(axonTree(:,1,1),axonTree(:,2,1),axonTree(:,3,1),'color',[1 0 0],'lineWidth',1);
% % plot the jittered axons
% for i = 1:Shuffles
% plot3(axonTreeJitter(:,1,i),axonTreeJitter(:,2,i),axonTreeJitter(:,3,i),'color',[0.8, 0.8, 0],'lineWidth',1);
% end
% 
% box on;
% axis([ 20000 140000 60000 250000 -60000 0]);
% plot( [20000, 40000], [70000, 70000],'-k' ) % insert 20um sclaebar
% daspect([1 1 1]); % make aspect ratio [1 1 1]


% calculate euclidean distance between axTree and denTree

clear Dist;
distThreshold = 1000; % set threshold for distance
%Dist = zeros(size(A,1),size(A,1),Shuffles);
for i = 1:1:Shuffles
Dist = pdist2(axonTreeJitter(:,:,i), denTree);
[temp1,temp2] = find(Dist>0 & Dist<distThreshold);
r(1:size(temp1,1),i) = temp1;
c(1:size(temp2,1),i) = temp2;
clear temp1;
clear temp2;
sprintf('Jitter iteration: %d',i)
end


% unique locations on axon

for i = 1:1:Shuffles
    if ~isempty(r(:,i))
        sprintf('Number of locations on the axon:%d \n Number of locations on dendrites: %d', length(r), length(c));
        sprintf('Number of unique locations on the axon:%d \n Number of unique locations on dendrites: %d', length(unique(r)), length(unique(c)));
        temp3 = find(diff(unique(r(:,i)))>1);
        diffUniq(1:size(temp3,1),i) = temp3;
        if ~isempty(diffUniq)
            uniquec(1:size( unique(c(:,i)),1),i) = unique(c(:,i));
            uniquer(1:size(unique(r(:,i)),1),i) = unique(r(:,i));
            %         s1 = scatter3(axonTree(uniquer(diffUniq),1),axonTree(uniquer(diffUniq),2),axonTree(uniquer(diffUniq),3),'Marker','o','MarkerFaceColor',[0 0 1],'MarkerEdgeColor','k');
            %         drawnow;
            %         s1Markers = s1.MarkerHandle;
            %         s1Markers.FaceColorData = uint8(255*[0;0;1;0.5]);  % Alpha=0.3 => 70% transparent red
            %         s1.SizeData = 50;
        end
        %unqSites(:,1:3,i) = [axonTreeJitter(uniquer(diffUniq(:,i)),1),axonTree(uniquer(diffUniq),2),axonTree(uniquer(diffUniq),3)];
    else
        %sprintf('no unique sites')
        uniquec(1,i) = 0;
        uniquer(1,i) = 0;
    end
    i
end

%unqSites = [axonTree(uniquer(diffUniq),1),axonTree(uniquer(diffUniq),2),axonTree(uniquer(diffUniq),3)];

% set figure to desired orientation

set (gca,'XTick',[], 'YTick',[],'ZTick', [],'Ydir','reverse');
view([-180,90]);
axis vis3d;
text(0.5,0.98,cellIDs{AxonTree},'Units','normalized');