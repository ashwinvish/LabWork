function q=getPoints(fileName,validNodeTypes)
if nargin<2; validNodeTypes=[-1:5]; end;
[nodes,edges,radii,nodeTypes,abort] = newReadSWCfile(fileName,validNodeTypes);

% extract positions
y = [nodes(edges(:,1),2) nodes(edges(:,2),2)]';
x = [nodes(edges(:,1),1) nodes(edges(:,2),1)]';
z = [nodes(edges(:,1),3) nodes(edges(:,2),3)]';
q = [x(:) y(:) -z(:)];
N = size(q,1);
% original scaling is pixels and planes 
%q = q.*[ones(N,2)*(0.005) ones(N,1)*(-0.045)]; % microns
%q = q.*[ones(N,2)*(0.005/0.045) -ones(N,1)]; % planes



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

 
