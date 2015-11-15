function [axonTree, axonTreeJitter] = UniqueSites(Axon, cellIDs, AxonTree, Display, jitter, jitterRadius)
%UNIQUESITES give the number of unique sites on the axon of a tree that are
%withing a certain radius of the other trees
%   Axon is a cell with node representation of all trees
%   cellIDs is a cell with IDs as strings for all cells in the datset
%   AxonTree is a scalar of the trees whose axon is being considered
%   Display ture or false to display the tree
%   jitter true or false to jitter or not in a random orientation
%   jitterRadius is the radius wihthin which to jitter the tree.

% load all axonal nodes
load CellAxons.mat

% Presynaptic Cell

AxnTree = Axon;
inducingNodes = eval([cellIDs{AxonTree},'_axon']);

% Plot the presynaptic cell ( not working )
if Display == true
    DisplayTree(AxnTree,[1],false,inducingNodes,[1 0.5 0.3]);
    hold on;
end

% consider only axon of the tree
axonTree = [];
jitterXYZ = [rand*jitterRadius, rand*jitterRadius, -1*rand*jitterRadius];
axonTreeJitter = [];

for kk=1:numel(AxnTree)
    children = AxnTree{kk}{2};
    for mm = 1:numel(children)
        if ismember(kk,inducingNodes)
            Axtempx=[AxnTree{children(mm)}{3}(1); AxnTree{children(mm)}{4}{1}(:,1); AxnTree{kk}{3}(1)];
            Axtempy=[AxnTree{children(mm)}{3}(2); AxnTree{children(mm)}{4}{1}(:,2); AxnTree{kk}{3}(2)];
            Axtempz=-[AxnTree{children(mm)}{3}(3); AxnTree{children(mm)}{4}{1}(:,3); AxnTree{kk}{3}(3)];
            if Display == true
                h1 = plot3(Axtempx,Axtempy,Axtempz,'color',[1 0 0],'lineWidth',2);
                hold on;
            end
            axonTree = [axonTree; Axtempx Axtempy Axtempz];
            
            if jitter == true
                AxtempxJitter = Axtempx+jitterXYZ(1);
                AxtempyJitter = Axtempy+jitterXYZ(2);
                AxtempzJitter = Axtempz+jitterXYZ(3);
                if Display == true
                    h1Jit = plot3(AxtempxJitter,AxtempyJitter,AxtempzJitter,'color',[0 1 0],'lineWidth',1);
                end
                axonTreeJitter = [axonTreeJitter; AxtempxJitter AxtempyJitter AxtempzJitter];
            end
        else
            continue;
        end
    end
end

if Display == true
    box on;
    axis([ 20000 140000 60000 250000 -60000 0]);
    plot( [20000, 40000], [70000, 70000],'-k' )     % insert 20um sclaebar
    daspect([1 1 1]);                               % make aspect ratio [1 1 1]
    set (gca,'XTick',[], 'YTick',[],'ZTick', [], 'Ydir','reverse');
    view([-180,90]);                             % xy view
end
end
