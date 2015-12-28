function [ AxonNodes ] = AxonQueryNodes(tree,queryNodes)
%AXONQUERYNODES returns the coordinates of the queryNodes
%   tree is the trees axon that is being queried
%	queryNodes are the nodes that are being queried
index = 1;
for i = 1:numel(tree)
    if ismember(i,queryNodes)
        AxonNodes(index,:) = tree{i}{3};
        index = index+1;
    end
end
end

