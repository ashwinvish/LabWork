function [tree] = treeVis(swcPath, displyFigure)
% TREEVIS is to visualize a tree, given an .swc file
%		swcPath is the path to the .swc file
%		displayFigure, ture or false

[nodes,edges,radii,nodeTypes,abort] = newReadSWCfile(swcPath,[-1 0 1 2 3 4 5]);
%resolution = [0.36,0.36,2]; LM resolution
resolution = [5,5,45];

nodes(:,1) = nodes(:,1) * resolution(1); % xresolution in nm
nodes(:,2) = nodes(:,2) * resolution(2); % yresolution in nm
nodes(:,3) = nodes(:,3) * resolution(3); % zresolition in nm


tree = generateIrreducibleDoubleLinkedTree(nodes,edges,pi*radii.^2);
validNodes = [1:numel(tree)]; colorString = 'red'; newFigure = true; inducingNodes = []; highlightedNodes = [];
if displayFigure == true
figure;hold;
for kk=1:numel(tree)
    children = tree{kk}{2};
    for mm=1:numel(children)
        tempx=[tree{children(mm)}{3}(1); tree{children(mm)}{4}{1}(:,1); tree{kk}{3}(1)]; %/0.397;
        tempy=[tree{children(mm)}{3}(2); tree{children(mm)}{4}{1}(:,2); tree{kk}{3}(2)]; %/0.397;
        tempz=-[tree{children(mm)}{3}(3); tree{children(mm)}{4}{1}(:,3); tree{kk}{3}(3)]; %/0.5;
        if ismember(kk,inducingNodes) && ismember(children(mm),inducingNodes)
            plot3(tempy,tempx,tempz,colorString);
        else
            if ismember(kk,validNodes)
                plot3(tempy,tempx,tempz,'black');%plot3(tempy,tempx,tempz,'white');
            end
        end
    end
end
for kk=1:numel(highlightedNodes)
    plot3(tree{highlightedNodes(kk)}{3}(2),tree{highlightedNodes(kk)}{3}(1),tree{highlightedNodes(kk)}{3}(3),'Marker','o','MarkerSize',5,'Color','r')
    %plot3(tree{highlightedNodes(kk)}{3}(2)/0.397,tree{highlightedNodes(kk)}{3}(1)/0.397,-tree{highlightedNodes(kk)}{3}(3)/0.5,'Marker','o','MarkerSize',5,'Color','r')
end
end

function tree = generateIrreducibleDoubleLinkedTree(nodes,edges,areas)
% remove the root nodes until the first node
[nodes,edges,areas] = removeImproperRootNodes(nodes,edges,areas);
% remove zero branches
[nodes,edges,areas] = removeZeroLengthEdges(nodes,edges,areas);
h = hist(edges(:,2),[1:size(nodes,1)]); irreducibleNodes = union(find(h~=1),1); % 1 is the root node
% put all the irreducible nodes at the beginning
for kk = 1:numel(irreducibleNodes)
    % swap on edges
    edges(edges==kk) = 0;
    edges(edges==irreducibleNodes(kk)) = kk;
    edges(edges==0) = irreducibleNodes(kk);
    % swap on nodes
    tmp = nodes(kk,:);
    nodes(kk,:) = nodes(irreducibleNodes(kk),:);
    nodes(irreducibleNodes(kk),:) = tmp;
    % swap on areas
    tmpA = areas(kk);
    areas(kk) = areas(irreducibleNodes(kk));
    areas(irreducibleNodes(kk)) = tmpA;
end
numelNodes = numel(irreducibleNodes);
%initialize tree with root as 1
tree{1}{1} = []; tree{1}{3} = nodes(1,:);
for kk = 1:numelNodes
    tree{kk}{2} = [];
end
for kk = 2:numelNodes
    path = nodes(kk,:);
    tmpParent = edges(find(edges(:,1)==kk),2);
    while tmpParent > numelNodes
        path = [path; nodes(tmpParent,:)];
        tmpParent = edges(find(edges(:,1)==tmpParent),2);
    end
    tree{kk}{1} = tmpParent; tree{tmpParent}{2} = [tree{tmpParent}{2} kk];
    tree{kk}{3} = nodes(kk,:);
    tree{kk}{4}{1} = [path; nodes(tmpParent,:)];
end

function [nodes,edges,areas] = removeImproperRootNodes(nodes,edges,areas)
root = 1; % node 1 is assigned as the root node
expandedEdges = [edges edges(:,1)];
child = expandedEdges([zeros(size(edges,1),1) edges==root]>0);
while numel(child) < 2
    nodes = [nodes(1:root-1,:); nodes(root+1:end,:)];
    areas = [areas(1:root-1); areas(root+1:end)];
    edges(any(edges==root,2),:) = [];
    edges(edges>root) = edges(edges>root) - 1;
    expandedEdges = [edges edges(:,1)];
    if child > root
        child = child-1;
    end
    root = child;
    child = expandedEdges([zeros(size(edges,1),1) edges==root]>0);
end

function [nodes,edges,areas] = removeZeroLengthEdges(nodes,edges,areas)
kk = 1;
while kk < size(edges,1)
    node1 = edges(kk,1); node2 = edges(kk,2);
    if norm(nodes(node1,:)-nodes(node2,:)) < 1e-10 % careful
        edges(edges==node1) = node2;
        edges(edges>node1) = edges(edges>node1)-1;
        edges(kk,:) = []; % remove edge
        nodes(node1,:) = []; % remove node
    else
        kk = kk + 1;
    end
end

function [nodes,edges,radii,nodeTypes,abort] = newReadSWCfile(fileName,validNodeTypes)
abort = false; nodes = []; edges = []; nodeTypes = [];
if nargin < 2
    validNodeTypes = [0 2 3 4];
end
validNodeTypes = setdiff(validNodeTypes,1);
[nodeID, nodeType, xPos, yPos, zPos, radii, parentNodeID] = textread(fileName, '%u%d%f%f%f%f%d','commentstyle', 'shell');
if ~any(parentNodeID==-1)
    disp(strcat('root not found in ',fileName));
    nodes = []; edges = []; radii = []; nodeTypes = []; abort = true;
    return;
end
nodeType(find(parentNodeID==-1))=1; % Every tree should start from a node of type 1 (soma)
firstSomaNode = find(nodeType == 1 & parentNodeID == -1, 1);
somaNodes = find(nodeType == 1); somaX = mean(xPos(somaNodes)); somaY = mean(yPos(somaNodes)); somaZ = mean(zPos(somaNodes)); somaRadius = mean(radii(somaNodes));
xPos(firstSomaNode) = somaX; yPos(firstSomaNode) = somaY; zPos(firstSomaNode) = somaZ; radii(firstSomaNode) = somaRadius;
parentNodeID(ismember(parentNodeID,somaNodes)) = firstSomaNode; % assign a single soma parent
nodesToDelete = setdiff(somaNodes,firstSomaNode); % delete all the soma nodes except for the firstSomaNode
nodeID(nodesToDelete)=[]; nodeType(nodesToDelete)=[]; xPos(nodesToDelete)=[]; yPos(nodesToDelete)=[]; zPos(nodesToDelete)=[]; radii(nodesToDelete)=[]; parentNodeID(nodesToDelete)=[];
for kk = 1:numel(nodeID)
    while ~any(nodeID==kk)
        nodeID(nodeID>kk) = nodeID(nodeID>kk)-1;
        parentNodeID(parentNodeID>kk) = parentNodeID(parentNodeID>kk)-1;
    end
end
validNodes = nodeID(ismember(nodeType,validNodeTypes));
additionalValidNodes = [];
for kk = 1:numel(validNodes)
    thisParentNodeID = parentNodeID(validNodes(kk)); thisParentNodeType = nodeType(thisParentNodeID);
    while ~ismember(thisParentNodeType,validNodeTypes)
        if thisParentNodeType == 1
            break;
        end
        additionalValidNodes = union(additionalValidNodes, thisParentNodeID); nodeType(thisParentNodeID) = validNodeTypes(1);
        thisParentNodeID = parentNodeID(thisParentNodeID); thisParentNodeType = nodeType(thisParentNodeID);
    end
end
validNodes = [firstSomaNode; validNodes; additionalValidNodes']; validNodes = unique(validNodes);
nodeID = nodeID(validNodes); nodeType = nodeType(validNodes); parentNodeID = parentNodeID(validNodes);
xPos = xPos(validNodes); yPos = yPos(validNodes); zPos = zPos(validNodes); radii = radii(validNodes);
for kk = 1:numel(nodeID)
    while ~any(nodeID==kk)
        nodeID(nodeID>kk) = nodeID(nodeID>kk)-1;
        parentNodeID(parentNodeID>kk) = parentNodeID(parentNodeID>kk)-1;
    end
end
nodes = [xPos yPos zPos];
edges = [nodeID parentNodeID];
edges(any(edges==-1,2),:) = [];
nodeTypes = unique(nodeType)';
