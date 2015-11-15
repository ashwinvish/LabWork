function [uniqueTree] = PotentialSites( axonalTree , dendriticTree,  AxonTree, DendTree)
%POTENTIALSITES Used to calculate potential sites between two trees
%   axonalTree is the axonal tree
%   dendriticTree is the dendritic tree
%   AxonTree is the cell number associated with the axonTree
%   DendTree is the cell number associated with the dendriticTree

if AxonTree == DendTree;
    uniqueTree = [];
    return;
end

clear Dist;
distThreshold = 1000;                                       % set threshold for distance
Dist = pdist2(axonalTree, dendriticTree);
[r,c] = find(Dist>0 & Dist<distThreshold);

% unique locations on axon

if ~isempty(r)
    diffUniq = find(diff(unique(r))>1);
    if ~isempty(diffUniq)
        uniquer = unique(r);
        uniqueTree = [axonalTree(uniquer(diffUniq),1),axonalTree(uniquer(diffUniq),2), axonalTree(uniquer(diffUniq),3)];
    else
        uniquer = 0;
        uniqueTree = [];
    end
else
    uniqueTree = [];
end
end




