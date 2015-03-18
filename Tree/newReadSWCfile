
function [nodes,edges,radii,nodeType,abort] = newReadSWCfile(fileName,validNodeTypes)
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
%nodeTypes = unique(nodeType)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nodes,edges,radii] = removeImproperRootNodes(nodes,edges,radii)
if size(nodes,1) < 2
  return;
end
% node 1 is assigned as the root node
root = 1;
expandedEdges = [edges edges(:,1)];
child = expandedEdges([zeros(size(edges,1),1) edges==root]>0);
while numel(child) < 2
  nodes = [nodes(1:root-1,:); nodes(root+1:end,:)];
  radii = [radii(1:root-1); radii(root+1:end)];
  edges(any(edges==root,2),:) = [];
  edges(edges>root) = edges(edges>root) - 1;
  expandedEdges = [edges edges(:,1)];
  if child > root
    child = child-1;
  end
  root = child;
  child = expandedEdges([zeros(size(edges,1),1) edges==root]>0);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

