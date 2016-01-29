function [Pot] = Shuffle(allTrees, cellIDs, cellIDsAlx, Display, jitterRadius, Shuffles )
%SHUFFLE shuffles a trees axon with respect to all other dendrites from
%remaining trees in the datset
%   allTrees is a cell with all the tress in the datset
%   cellIDS is a cell containign IDs of all cells as strings
%   cellIDsAlx is the Ids of a unique population from CellIDs
%   Display true or false
%   jitterRadius is radius by which a tree should be translated in a random
%   orinetation
%   Shuffles is the number of times the tree should be moved.

load CellAxons.mat

for i = 1:numel(cellIDs)
    
    if ismember(cellIDs(i), cellIDsAlx) == 1
        
        AxonTree = i
        thisAxonTree = allTrees{AxonTree};
        % pre allocate
        [A] =  UniqueSites(thisAxonTree,cellIDs,AxonTree,false,true, jitterRadius);
        % axonTree is (m,n,p) sized. m is the number of nodes in the tree,
        % n is 3 coordinates, p is the same size of number of shuffles
        axonTree = zeros(size(A,1),3);
        axonTreeJitter = zeros(size(A,1),3,Shuffles);
        
        % polulate jittered axonal tree
        for j = 1:1:Shuffles;
            [axonTree(1:size(A,1),1:3), axonTreeJitter(1:size(A,1),1:3,j)]= UniqueSites(thisAxonTree,cellIDs,AxonTree,false,true, jitterRadius);
        end
        %  if the cell does not have an axon reuturn
        if isempty(axonTree)
            diffUniq = [];
            %return
            continue;
        end
        
        % consider dendritic trees of all remaining cells
        
        for kk = 1:Shuffles
            iteration = kk;
            uniqSites = [];
            
            for ii = 1:numel(cellIDs)
                denTree = [];
                DenTree = allTrees{ii};                                         % iterating through all trees
                validNodes =  eval([cellIDs{ii},'_axon']);                       % keeping track of axonal nodes
                if isempty(validNodes)                                           % if no axon check
                    validNodes = 1:numel(DenTree);
                end                                                              % if no axon then iterate through all nodes
                for jj = 1:numel(DenTree)
                    children = DenTree{jj}{2};                                   % Consider childern of jj node
                    for nn = 1:numel(children)                                   % iterate over all children
                        DnTempx = [DenTree{children(nn)}{3}(1); DenTree{children(nn)}{4}{1}(:,1)];
                        DnTempy = [DenTree{children(nn)}{3}(2); DenTree{children(nn)}{4}{1}(:,2)];
                        DnTempz = -[DenTree{children(nn)}{3}(3); DenTree{children(nn)}{4}{1}(:,3)];
                        if length(validNodes)< length(DenTree) && ismember(jj,validNodes) % dont plot the axon
                            continue;
                        else
                            if Display == true
                                h2 = plot3(DnTempx,DnTempy,DnTempz,'-','color',[0.8 0.8 0.8]);
                                %h2.Color(4) = 0.2;
                            end
                        end
                        denTree = [denTree; DnTempx DnTempy DnTempz];            % populate with all dendritic trees of all trees
                    end
                end
                
                [TempuniqSites] = PotentialSites(axonTreeJitter(:,:,kk), denTree, i,ii);
                [uniqSites] = [uniqSites;TempuniqSites];
            end
            Pot(i,kk) = size(uniqSites,1);
        end
        
    end
end
