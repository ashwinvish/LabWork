function [ tree,rawLength ] = generateTree( nodes,edges,radii, uniformCSAreas )
% generateTree displays a tree given the graph
%   Detailed explanation goes here

%[nodes,edges,radii] = removeZeroLengthEdges(nodes,edges,radii); % remove zero branches
areas = pi*radii.^2;
% 1 is assumed to be the root node. edges is nx2: each row is (child,
% parent)
h = hist(edges(:,2),[1:size(nodes,1)]); 
irreducibleNodes = union(find(h~=1),1); % 1 is the root node

% put all the irreducible nodes at the beginning
for kk = 1:numel(irreducibleNodes)
    % swap on edges
    edges(edges==kk) = 0; edges(edges==irreducibleNodes(kk)) = kk; edges(edges==0) = irreducibleNodes(kk);
    % swap on nodes
    tmp = nodes(kk,:); nodes(kk,:) = nodes(irreducibleNodes(kk),:); nodes(irreducibleNodes(kk),:) = tmp;
    % swap on areas
    tmpA = areas(kk); areas(kk) = areas(irreducibleNodes(kk)); areas(irreducibleNodes(kk)) = tmpA;
end

numelNodes = numel(irreducibleNodes);

%initialize tree with root as 1
tree{1}{1} = [];
tree{1}{3} = nodes(1,:); 
tree{1}{4}{1} = nodes(1,:);
tree{1}{4}{2} = 0; 
tree{1}{4}{3} = [[0;0;0] [0;0;0]];
tree{1}{4}{4} = 0;

for kk = 1:numelNodes
    tree{kk}{2} = [];
end

rawLength = 0;

for kk = 2:numelNodes
    tmpParent = edges(find(edges(:,1)==kk),2); % assume that the edges are ordered pairs: (child, parent)
    path = nodes(kk,:);
    if uniformCSAreas
        pathAreas = 1;
    else
    
        pathAreas = (areas(kk)+areas(tmpParent)+sqrt(areas(kk)*areas(tmpParent)))/3;
    end
    
    while tmpParent > numelNodes
        newTmpParent = edges(find(edges(:,1)==tmpParent),2);
        path = [path; nodes(tmpParent,:)];
       
            pathAreas = [pathAreas; (areas(tmpParent)+areas(newTmpParent)+sqrt(areas(tmpParent)*areas(newTmpParent)))/3]; % now modeled as a cylinder
        tmpParent = newTmpParent;
    end
    
    path = [path; nodes(tmpParent,:)];
    rawPathLengths = sqrt(sum(diff(path,1,1).^2,2));
    rawLength = rawLength + sum(rawPathLengths);
    
    tree{kk}{1} = tmpParent;
    tree{tmpParent}{2} = [tree{tmpParent}{2} kk];
    tree{kk}{3} = nodes(kk,:);
    tree{kk}{4}{1} = path;
    tree{kk}{4}{2} = rawPathLengths;
    tree{kk}{4}{4} = pathAreas;
end

function [nodes,edges,radii] = removeZeroLengthEdges(nodes,edges,radii)
kk = 1;
while kk < size(edges,1)
    node1 = edges(kk,1); node2 = edges(kk,2);
    if norm(nodes(node1,:)-nodes(node2,:)) < 1e-10 % BE CAREFUL / ARBITRARY
        edges(edges==node1) = node2;
        edges(edges>node1) = edges(edges>node1)-1;
        edges(kk,:) = []; % remove edge
        nodes(node1,:) = []; % remove node
        radii(node1) = []; % remove area
    else
        kk = kk + 1;
    end
end



