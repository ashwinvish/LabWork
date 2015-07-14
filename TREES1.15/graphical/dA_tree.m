% DA_TREE   Plots the adjacency matrix of a tree.
% (trees package)
% 
% HP = dA_tree (intree, color, DD, xyscale, options)
% ---------------------------------------------------------
%
% displays the adjacency matrix of a tree by filling in an N x N square with
% 1s and 0s if N<50 and by black dots if tree is bigger
%
% Input
% -----
% - intree::integer:index of tree in trees structure or structured tree
%                     or simply an adjacency matrix
% - color::RGB 3-tupel: RGB values {DEFAULT black [0 0 0]}
% - DD:: XY-tupel or XYZ-tupel: coordinates offset {DEFAULT no offset [0,0,0]}
% - xyscale::value: scaling factor for X and Y  {DEFAULT no scaling ==1}
% - options::string: {DEFAULT ''}
%     '-t' : force text output
%     '-g' : force graphical output
%
% Output
% ------
% - HP::handles: depending on options HP links to the graphical objects.
%
% Example
% -------
% dA_tree (sample2_tree, [1 0 0]); axis off;
% dA_tree (sample_tree)
%
% See also start_trees
% Uses X,Y,Z
%
% the TREES toolbox: edit, visualize and analyze neuronal trees
% Copyright (C) 2009  Hermann Cuntz

function HP = dA_tree (intree, color, DD, xyscale, options)

% trees : contains the tree structures in the trees package
global trees

if (nargin < 1)||isempty(intree),
    intree = length(trees); % {DEFAULT tree: last tree in trees cell array}
end;

% use only directed adjacency for this function (obviously..)
if ~isstruct (intree),
    if numel (intree) > 1,
        dA = intree;
    else
        ver_tree (intree); % verify that input is a tree structure
        dA = trees {intree}.dA;
    end
else
    ver_tree (intree); % verify that input is a tree structure
    dA = intree.dA;
end

if (nargin<2)||isempty(color),
    color = [0 0 0]; % {DEFAULT: black}
end

if (nargin<3)||isempty(DD),
    DD = [0 0 0]; % {DEFAULT 3-tupel: no spatial displacement from the root}
end
if length(DD)<3,
    DD = [DD zeros(1, 3 - length (DD))]; % append 3-tupel with zeros
end

if (nargin<4)||isempty(xyscale),
    xyscale = 1; % {DEFAULT: no scaling}
end

if (nargin <5)||isempty(options),
    options = ''; % {DEFAULT: no option}
end

N = size (dA, 1); % number of nodes in tree
hold on;
[Y, X] = ind2sub (size (dA), find (dA)); % X and Y positions of non-zero entries of dA
X = xyscale * X + DD (1);
Y = xyscale * Y - DD (2);
if ((N < 50) && isempty(strfind (options, '-g'))) || ~isempty(strfind (options, '-t')),
    HP1 = line (repmat ( xyscale * (0.5 : 1 : N + 0.5), 2, 1)     + DD (1), ...
        repmat         (-xyscale * [0.5 N+0.5]',        1, N + 1) + DD (2));
    set (HP1, 'color', color);
    HP2 = line (repmat ( xyscale * [0.5 N+0.5]', 1, N + 1)    + DD (1), ...
        repmat         (-xyscale * (0.5 : 1 : N + 0.5), 2, 1) + DD (2));
    set (HP2, 'color', color);
    HT1 = text (xyscale * (1 : N) + DD (1), zeros (1, N)  + DD (2), num2str ((1 : N)'));
    set (HT1, 'horizontalalignment', 'center', 'color', color);
    HT2 = text (zeros (1, N) + DD (1), -xyscale * (1 : N) + DD (2), num2str ((1 : N)'));
    set (HT2, 'verticalalignment', 'middle', 'color', color);
    HT3 = text (X, -Y, '1');
    set (HT3,'horizontalalignment', 'center', 'verticalalignment', 'middle', 'color', color);
    HP3 = patch ((repmat (X, 1, 5) + repmat (xyscale * [-.5 .5 .5 -.5 -.5], length (X), 1))', ...
        (repmat         (-Y, 1, 5) + repmat (xyscale * [.5 .5 -.5 -.5 .5],  length (X), 1))', color);
    set (HP3, 'facealpha', 0.2, 'linestyle', 'none');
    [Y, X] = ind2sub (size (dA), find (~dA)); % the opposite... fill with "0"
    X = xyscale * X + DD (1); Y = xyscale * Y - DD (2);
    HT4 = text (X, -Y, '0');
    set (HT4,'horizontalalignment', 'center', 'verticalalignment', 'middle', 'color', color);
    HP  = [HP1; HP2; HP3; HT1; HT2; HT3; HT4];
else
    HP1 = line (repmat (xyscale * [0.5 N+0.5]  + DD (1), 2, 1), ...
        repmat        (-xyscale * [0.5 N+0.5]' + DD (2), 1, 2));
    set (HP1, 'color', color);
    HP2 = line (repmat (xyscale * [0.5 N+0.5]' + DD (1), 1, 2),...
        repmat        (-xyscale * [0.5 N+0.5]  + DD (2), 2, 1));
    set(HP2, 'color', color);
    HT1 = text (xyscale * [1 N]+ DD (1), zeros (1, 2)+ DD (2), num2str ([1 N]'));
    set (HT1, 'verticalalignment', 'bottom', 'horizontalalignment', 'center', 'color', color);
    HT2 = text (zeros (1, 2) + DD (1), -xyscale * [1 N]+ DD (2), num2str ([1 N]'));
    set (HT2, 'horizontalalignment', 'right', 'color', color);
    HP3 = patch ((repmat (X, 1, 5) + repmat (xyscale * [-.5 .5 .5 -.5 -.5], length (X), 1))',...
        (repmat         (-Y, 1, 5) + repmat (xyscale * [.5 .5 -.5 -.5 .5],  length (X), 1))', color);
    set (HP3, 'facealpha', 0.2, 'linestyle', 'none');
    HP4 = plot (X, -Y, 'ko');
    set (HP4, 'color', color, 'Markerfacecolor', color, 'markersize', 4);
    HP  = [HP1; HP2; HP3; HP4; HT1; HT2];
end
axis equal
