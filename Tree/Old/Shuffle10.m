% shuffle axon tree

for i = 1:numel(cellIDs)
    if ismember(cellIDs(i), cellIDsAlx) == 1
        AxonTree = i
        Shuffles = 1;
        
        
        [A,B] =  UniqueSites10(allTrees,cellIDs,AxonTree,false,true);
        axonTree = zeros(size(A,1),3,Shuffles);
        axonTreeJitter = zeros(size(A,1),3,Shuffles);
        
        for j = 1:1:Shuffles;
            [axonTree(1:size(A,1),1:3,j), axonTreeJitter(1:size(A,1),1:3,j)]= UniqueSites10(allTrees,cellIDs,AxonTree,false,true, 10000);
        end
        
        if isempty(axonTree)
            diffUniq = [];
            return
        end
        
        
        % consider dendritic trees of all remaining cells
        
        remTrees = 1:numel(cellIDs);
        remTrees(AxonTree) = [];
        denTree = [];
        %        hold on;
        
        for ii = 1:numel(remTrees)
            DenTree = allTrees{remTrees(ii)};                % iterating through all trees
            validNodes =  eval([cellIDs{ii},'_axon']);       % keeping track of axonal nodes
            if isempty(validNodes)                           % if no aoxon check
                validNodes = 1:numel(DenTree);
            end                                             % if no axon then iterate through all nodes
            for jj = 1:numel(DenTree)
                children = DenTree{jj}{2};                   % Consider childern of jj node
                for nn = 1:numel(children)                   % iterate over all children
                    DnTempx = [DenTree{children(nn)}{3}(1); DenTree{children(nn)}{4}{1}(:,1)];
                    DnTempy = [DenTree{children(nn)}{3}(2); DenTree{children(nn)}{4}{1}(:,2)];
                    DnTempz = -[DenTree{children(nn)}{3}(3); DenTree{children(nn)}{4}{1}(:,3)];
                    if length(validNodes)< length(DenTree) && ismember(jj,validNodes)
                        continue;
                    else
                        %h2 = plot3(DnTempx,DnTempy,DnTempz,'-','color',[0.8 0.8 0.8]);
                        %h2.Color(4) = 0.2;
                    end
                    denTree = [denTree; DnTempx DnTempy DnTempz];
                end
            end
        end
        
        
        % calculate euclidean distance between axTree and denTree
        
        clear Dist;
        distThreshold = 1000; % set threshold for distance
        Dist = zeros(size(axonTreeJitter(:,:,1),1), size(denTree,1));
        
        for kk = 1:1:Shuffles
            Dist = pdist2(axonTreeJitter(:,:,kk), denTree);
            [temp1,temp2] = find(Dist>0 & Dist<distThreshold);
            r(1:size(temp1,1),kk) = temp1;
            c(1:size(temp2,1),kk) = temp2;
            clear temp1;
            clear temp2;
            sprintf('Jitter iteration: %d',kk)
        end
        
        
        % find unique locations on axon
        clear diffUniq;
        clear temp3;
        clear uniquer;
        clear uniquec;
        
        for ll = 1:1:Shuffles
            if ~isempty(r(:,ll))
                sprintf('Number of locations on the axon:%d \n Number of locations on dendrites: %d', length(r), length(c));
                sprintf('Number of unique locations on the axon:%d \n Number of unique locations on dendrites: %d', length(unique(r)), length(unique(c)));
                temp3 = find(diff(unique(r(:,ll)))>1);
                diffUniq(1:size(temp3,1),ll) = temp3;
            else
                diffUniq(1:size(temp3,1),ll) = [];
            end
            ll
        end
        %str =  sprintf('JitteredPotSynapse_%s_20um',cell2mat(cellIDs(AxonTree)));
%         for mm = 1:1:Shuffles
%             Pot(AxonTree,mm) = size(diffUniq((diffUniq(:,mm)~=0)),1);
%         end
       % save('10umJitter.mat', 'Pot');
    end
end
