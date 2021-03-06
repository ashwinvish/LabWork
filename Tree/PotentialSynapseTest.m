
% Plot axon of presynaptic cells and dendrites of all postsynaptic cells
% along with potential synapses

load CellAxons.mat
jitter = true;                                                               % ture or false
jitterRadius = 5000;                                                         % radius in nm

figure();
hold on;

% Presynaptic Cell
A = [4]; 
AxnTree = allTrees{A};

 if jitter == true
     [axonTree,axonTreeJitter] = UniqueSites(AxnTree,cellIDs,A,true,true,jitterRadius);
 else
     [axonTree,axonTreeJitter] = UniqueSites(AxnTree,cellIDs,A, true,true, (jitterRadius-jitterRadius));
 end

% consider dendritic trees of all remaining cells

uniq= [];
uniqJitter = [];


for ii = 1:numel(cellIDs)
    denTree = [];
    DenTree = allTrees{ii};                                        % iterating through all trees
    validNodes =  eval([cellIDs{ii},'_axon']);                     % keeping track of axonal nodes
    if isempty(validNodes)                                                   % if no axon check
        validNodes = 1:numel(DenTree);
    end                                                                      % if no axon then iterate through all nodes
    for jj = 1:numel(DenTree)
        children = DenTree{jj}{2};                                           % Consider childern of jj node
        for nn = 1:numel(children)                                           % iterate over all children
            DnTempx = [DenTree{children(nn)}{3}(1); DenTree{children(nn)}{4}{1}(:,1)];
            DnTempy = [DenTree{children(nn)}{3}(2); DenTree{children(nn)}{4}{1}(:,2)];
            DnTempz = -[DenTree{children(nn)}{3}(3); DenTree{children(nn)}{4}{1}(:,3)];
            if length(validNodes)< length(DenTree) && ismember(jj,validNodes)   
                continue;
            else
                h2 = plot3(DnTempx,DnTempy,DnTempz,'-','color',[0.8 0.8 0.8]);
            end
            denTree = [denTree; DnTempx DnTempy DnTempz];
        end
    end

    
    [Tempuniq] = PotentialSites(axonTree,denTree,A,ii);
    [uniq] = [uniq;Tempuniq];
    clear Tempuniq
    [TempuniqJitter] = PotentialSites(axonTreeJitter,denTree,A,ii);
    [uniqJitter] = [uniqJitter;TempuniqJitter];
    clear TempuniqJitter;
    
    
    clear DenTree;
    clear validNodes;
    clear children;
    clear denTree;
    
end

scatter3(uniq(:,1),uniq(:,2),uniq(:,3),'Marker','o','MarkerFaceColor','b'); % plot all potential synapses
scatter3(uniqJitter(:,1),uniqJitter(:,2), uniqJitter(:,3), 'Marker','o','MarkerFaceColor','g'); % plot all jittered synapses
str = sprintf('Potential synapses %d \n Jittered potential synapses %d', size(uniq,1),size(uniqJitter,1));
title(str);

% size(uniq,1)
% size(uniqJitter,1)

% calculate euclidean distance between axTree and denTree


