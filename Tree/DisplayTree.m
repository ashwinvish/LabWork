function DisplayTree( tree, visualNodes, newFig, axonHighlight, color, presynapses, postsynapses )
%DISPLAYTREE Displays the neuron tree
%   tree is the neuron
%   newFig plot in new figure, true or false
%   axonHighlight are the nodes that represent the axon
%   color is the 3 channel color by which the tree is represented
%   presynapses are the presynaptic locations
%   postsynapses are the postsynaptic locations

load CellAxons.mat;

if nargin < 7
    postsynapses = [];
    if nargin < 6
        presynapses = [];
        %if nargin == 4
        %   axonHighlightColor  = abs(color-0.3) ;
        if nargin < 5
            color = [rand, rand , rand];
            if nargin < 4
                axonHighlight = [];
                if nargin < 3
                    newFig = true;
                    if nargin < 2
                        visualNodes = [1];
                    end
                end
            end
        end
    end
end

axonHighlightColor  = abs(color-0.3) ;
treeVisualizer(tree, visualNodes, axonHighlight,[{postsynapses},{presynapses}],newFig,{color axonHighlightColor}, 1:numel(tree), false);
end

