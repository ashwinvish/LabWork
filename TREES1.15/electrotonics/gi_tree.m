% GI_TREE   Axial conductances of the segments of a tree.
% (trees package)
% 
% gi = gi_tree (intree, options)
% ------------------------------
%
% returns the axial conductances of all elements [in Siemens].
%
% Input
% -----
% - intree::integer:index of tree in trees or structured tree
% - options::string: {DEFAULT: ''}
%     '-s'  : show
%
% Output
% -------
% gi::Nx1 vector: axial conductance values of each segment
%
% Example
% -------
% gi_tree (sample_tree, '-s')
%
% See also gm_tree
% Uses cvol_tree ver_tree Ri
%
% the TREES toolbox: edit, visualize and analyze neuronal trees
% Copyright (C) 2009  Hermann Cuntz

function gi = gi_tree (intree, options)

% trees : contains the tree structures in the trees package
global trees

if (nargin < 1)||isempty(intree),
    intree = length (trees); % {DEFAULT tree: last tree in trees cell array}
end;

ver_tree (intree); % verify that input is a tree structure

% use only axial resistance vector/value for this function
if ~isstruct (intree),
    Ri = trees {intree}.Ri;
else
    Ri = intree.Ri;
end

if (nargin < 2)||isempty(options),
    options = ''; % {DEFAULT: no option}
end;

Hlov = 1 ./ (cvol_tree (intree) * 10000);
% conversion cvol from 1/um to 1/cm Hlov is in [cm]
gi   = Hlov ./ Ri;

if strfind (options, '-s'), % show option
    ipart = find (gi < 0.0099); % single out non-0-length segments
    clf; shine; hold on; plot_tree (intree, gi, [], ipart); colorbar;
    title  ('axial conductances [S]');
    xlabel ('x [\mum]'); ylabel ('y [\mum]'); zlabel ('z [\mum]');
    view (2); grid on; axis image;
end

