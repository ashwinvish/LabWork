function [axonTree, axonTreeJitter] = UniqueSites(AllTrees, cellIDs, AxonTree, DisplayTree, jitter)
%UNTITLED23 Summary of this function goes here
%   Detailed explanation goes here


load CellAxons.mat

% Plot axon of presynaptic cells and dendrites of all postsynaptic cells

% %AllTrees = AllTrees;
% AxonTree = [5];
% % Presynaptic Cell
%AxnTree = AllTrees{AxonTree};

A = [5]; 

AxnTree = AllTrees{A};

% Plot the presynaptic cell

if DisplayTree == true
    hold on;
    DisplayTree(AxnTree,[1],true,eval([cellIDs{AxonTree},'_axon']),[1 0.5 0.3]);
    %DisplayTree(AxnTree,[1],true,eval([cellIDs{AxonTree},'_axon']),[1 0.5 0.3]);
end


% consider only axon of the tree

inducingNodes = eval([cellIDs{AxonTree},'_axon']);
axonTree = [];
jitterRadius = 30*1000 ; % in nm
jitterXYZ = [rand*jitterRadius, rand*jitterRadius, -1*rand*jitterRadius];
axonTreeJitter = [];

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

if DisplayTree == true
    box on;
    axis([ 20000 140000 60000 250000 -60000 0]);
    plot( [20000, 40000], [70000, 70000],'-k' ) % insert 20um sclaebar
    daspect([1 1 1]); % make aspect ratio [1 1 1]
end


end
