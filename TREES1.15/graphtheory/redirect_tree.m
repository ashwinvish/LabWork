% REDIRECT_TREE   Set root to new point and redirect tree graph.
% (trees package)
% 
% [tree, order] = redirect_tree (intree, istart, options)
% -------------------------------------------------------
%
% changes the direction of the adjacency matrix so that arrows show away
% from element istart (which becomes the first elements)
% ! redirect only makes sense on terminal or continuation nodes !
% ! because otherwise a trifurcation occurs                   !
%
% Input
% -----
% - intree::integer: index of tree in trees or structured tree
% - istart::number: starting node {DEFAULT: last node}
% - options::string: {DEFAULT: ''}
%     '-s' : show (green - old index, red - new index)
%
% Outputs
% -------
% if no output is declared the tree is changed in trees
% - tree:: structured output tree
% - order::Nx1 vector: vector of new indices
%
% Example
% -------
% redirect_tree (sample2_tree, 5, '-s')
%
% See also sortBCT_tree sortLO_tree
% Uses PL_tree ver_tree dA
%
% the TREES toolbox: edit, visualize and analyze neuronal trees
% Copyright (C) 2009  Hermann Cuntz

function varargout = redirect_tree (intree, istart, options)

% trees : contains the tree structures in the trees package
global trees

if (nargin < 1)||isempty(intree),
    intree = length(trees); % {DEFAULT tree: last tree in trees cell array}
end

ver_tree(intree); % verify that input is a tree structure

% use full tree for this function
if ~isstruct(intree),
    tree = trees{intree};
else
    tree = intree;
end

if (nargin < 2)||isempty(istart),
    istart = size(tree.dA,1); % {DEFAULT index: last node in tree}
end
    
if (nargin < 3)||isempty(options),
    options = ''; % {DEFAULT: no option}
end

if sum(tree.dA(:,istart))==2,
    warning('TREES:BCTinconsistency','branching point! => trifurcation will occur');
end

% simple use of the adjacency matrix again
A = tree.dA + tree.dA';
ward = 1;
W = A (:,istart);
PL = W;
resW = W;
% maximum depth by 2x maximum path length for when root is in center:
maxPL = max(PL_tree (tree))*2;
while ((sum (sum (resW == 1)) ~= 0) && (ward <= maxPL)),
    ward = ward + 1;
    resW = A*resW; % use adjacency matrix to walk through tree
    PL = PL + ward.*(resW == 1);
end
PL (istart)=0;
[t ii] = sort (PL);

% change the trees-structure according to the new order:
tree.dA = tril(A(ii,ii));
% in all vectors of form Nx1
S = fieldnames (tree);
for ward = 1:length(S),
    if ~strcmp(S{ward},'dA'),
        vec = tree.(S{ward});
        if isvector(vec) && (numel(vec) == size (tree.dA, 1)),
            tree.(S{ward}) = tree.(S{ward})(ii);
        end
    end
end

if strfind(options,'-s'), % show option
    clf; hold on; shine; HP = plot_tree (intree); set(HP, 'facealpha',0.2);
    T = vtext_tree (intree, [], [0 1 0], [-2 3 5]); set (T, 'fontsize', 14);
    T = vtext_tree (tree, [], [], [0 0 5]); set (T, 'fontsize', 14);
    title ('redirect to other root node');
    HP(1) = plot(1, 1, 'g-'); HP(2) = plot(1, 1, 'r-');
    legend (HP, {'before', 'after'});
    set (HP, 'visible', 'off');
    xlabel ('x [\mum]'); ylabel ('y [\mum]'); zlabel ('z [\mum]');
    view (2); grid on; axis image;
end

if (nargout >0||(isstruct(intree)))
    varargout{1} = tree; % if output is defined then it becomes the tree
else
    trees{intree} = tree; % otherwise add to end of trees cell array
end

if (nargout>1)
    varargout{2} = ii;
end
