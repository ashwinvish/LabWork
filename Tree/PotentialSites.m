function [uniqueTree] = PotentialSites( allTrees , cellIDs,  AxonTree, DendriticTree, DisplayTree, jitter)
%UNTITLED26 Summary of this function goes here
%   Detailed explanation goes here


load CellAxons.mat


% Presynaptic Cell
AxnTree = allTrees{AxonTree};

% Plot the presynaptic cell


% consider only axon of the tree

inducingNodes = eval([cellIDs{AxonTree},'_axon']);
axonTree = [];

if jitter == true
    jitterRadius = 10*1000 ; % in nm
    jitterXYZ = [rand*jitterRadius, rand*jitterRadius, -1*rand*jitterRadius];
    axonTreeJitter = [];
end

for kk=1:numel(AxnTree)
    children = AxnTree{kk}{2};
    for mm = 1:numel(children)
        if ismember(kk,inducingNodes)
            Axtempx=[AxnTree{children(mm)}{3}(1); AxnTree{children(mm)}{4}{1}(:,1); AxnTree{kk}{3}(1)];
            Axtempy=[AxnTree{children(mm)}{3}(2); AxnTree{children(mm)}{4}{1}(:,2); AxnTree{kk}{3}(2)];
            Axtempz=-[AxnTree{children(mm)}{3}(3); AxnTree{children(mm)}{4}{1}(:,3); AxnTree{kk}{3}(3)];
            if DisplayTree == true
                h1 = plot3(Axtempx,Axtempy,Axtempz,'color',[1 0 0],'lineWidth',1);
                %h1.Color(4) = 0.2;                                             %  transparency 0-1
                hold on;
            end
            axonTree = [axonTree; Axtempx Axtempy Axtempz];
            if jitter == true
                AxtempxJitter = Axtempx+jitterXYZ(1);
                AxtempyJitter = Axtempy+jitterXYZ(2);
                AxtempzJitter = Axtempz+jitterXYZ(3);
                if DisplayTree == true
                    h1Jit = plot3(AxtempxJitter,AxtempyJitter,AxtempzJitter,'color',[0.5 0 0],'lineWidth',1);
                end
                
                axonTreeJitter = [axonTreeJitter; AxtempxJitter AxtempyJitter AxtempzJitter];
            end
            
        else
            continue;
        end
    end
end

if isempty(axonTree)
    uniqueTree = [];
    return;
end


denTree = [];

DenTree = allTrees{DendriticTree};                          % iterating through all trees
validNodes =  eval([cellIDs{DendriticTree},'_axon']);       % keeping track of axonal nodes
if isempty(validNodes)                                      % if no aoxon, iterate over dendritic nodes
    validNodes = 1:numel(DenTree);
end                                                          % if no axon then iterate through all nodes
for jj = 1:numel(DenTree)
    children = DenTree{jj}{2};                              % Consider childern of jj node
    for nn = 1:numel(children)                              % iterate over all children
        DnTempx = [DenTree{children(nn)}{3}(1); DenTree{children(nn)}{4}{1}(:,1)];  % DenTree{nn}{3}(1)];
        DnTempy = [DenTree{children(nn)}{3}(2); DenTree{children(nn)}{4}{1}(:,2)];  % DenTree{nn}{3}(2)];
        DnTempz = -[DenTree{children(nn)}{3}(3); DenTree{children(nn)}{4}{1}(:,3)]; % DenTree{nn}{3}(3)];
        if length(validNodes)< length(DenTree) && ismember(jj,validNodes)
            continue;
        else
            h2 = plot3(DnTempx,DnTempy,DnTempz,'-','color',[0.8 0.8 0.8]);%,'lineWidth',0.5);
            %h2.Color(4) = 0.2;
        end
        denTree = [denTree; DnTempx DnTempy DnTempz];
    end
end

box on;
axis([ 20000 140000 60000 250000 -60000 0]);
plot( [20000, 40000], [70000, 70000],'-k' )                 % insert 20um sclaebar
daspect([1 1 1]);                                           % make aspect ratio [1 1 1]
view(-180,90);
set (gca,'XTick',[], 'YTick',[],'ZTick', [],'Ydir','reverse');

% if AxonTree == DendriticTree;
%     uniqueTree = [];
%     return;
% end


clear Dist;
distThreshold = 1000;                                       % set threshold for distance
Dist = pdist2(axonTree, denTree);
[r,c] = find(Dist>0 & Dist<distThreshold);


% unique locations on axon

if ~isempty(r)
    sprintf('Number of locations on the axon:%d \n Number of locations on dendrites: %d', length(r), length(c));
    %diffUniq = find(diff(unique(c))>1);
    diffUniq = find(diff(unique(r))>1);
    if ~isempty(diffUniq)
        %sprintf('Number of unique locations on the axon:%d \n Number of unique locations on dendrites: %d', length(unique(r)), length(unique(c)))
        uniquec = unique(c);
        uniquer = unique(r);
        %s1 = scatter3(axonTree(uniquer(diffUniq),1),axonTree(uniquer(diffUniq),2),axonTree(uniquer(diffUniq),3),'Marker','o','MarkerFaceColor',[0 0 1],'MarkerEdgeColor','k');
        uniqueTree = [axonTree(uniquer(diffUniq),1),axonTree(uniquer(diffUniq),2), axonTree(uniquer(diffUniq),3)];
        %drawnow;
        %s1Markers = s1.MarkerHandle;
        %s1Markers.FaceColorData = uint8(255*[0;0;1;0.5]);  % Alpha=0.3 => 70% transparent red
        %s1.SizeData = 50;
        %title(sprintf('Potential Synapses %d', size(uniqueTree,1)));
        
    else
        %sprintf('no unique sites')
        uniquec = 0;
        uniquer = 0;
        uniqueTree = [];
    end
else
    uniqueTree = [];
end
end




