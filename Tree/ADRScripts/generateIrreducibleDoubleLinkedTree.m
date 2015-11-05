function [tree,rawLength,specialPositions] = generateIrreducibleDoubleLinkedTree(fileName,validNodeTypes,SpecialNodeType,uniformCSAreas)

[nodes,edges,radii,nodeType,temp2] = newReadSWCfile(fileName,validNodeTypes);
specialPositions = nodes(ismember(nodeType,SpecialNodeType),:);
areas = pi*radii.^2;
 % 1 is assumed to be the root node. edges is nx2: each row is [child parent]
 h = hist(edges(:,2),[1:size(nodes,1)]); irreducibleNodes = union(find(h~=1),1); % 1 is the root node
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
 tree{1}{1} = []; tree{1}{3} = nodes(1,:); tree{1}{4}{1} = nodes(1,:); tree{1}{4}{2} = 0; tree{1}{4}{3} = [[0;0;0] [0;0;0]]; tree{1}{4}{4} = 0; 
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
     if uniformCSAreas
       pathAreas = [pathAreas; 1];
     else
       pathAreas = [pathAreas; (areas(tmpParent)+areas(newTmpParent)+sqrt(areas(tmpParent)*areas(newTmpParent)))/3]; % now modeled as a cylinder
     end
     tmpParent = newTmpParent;
   end
   path = [path; nodes(tmpParent,:)]; rawPathLengths = sqrt(sum(diff(path,1,1).^2,2)); rawLength = rawLength + sum(rawPathLengths);
   tree{kk}{1} = tmpParent; tree{tmpParent}{2} = [tree{tmpParent}{2} kk];
   tree{kk}{3} = nodes(kk,:);
   tree{kk}{4}{1} = path; tree{kk}{4}{2} = rawPathLengths; tree{kk}{4}{4} = pathAreas;
 end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

