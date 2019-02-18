function [pathLength] = PathLengthToCoordinate(pointCloud,inTree)
%% Calculate pathlength from root node to A
% A is mx3 vector
% inTree is the tree structure associated with that A.

B = [inTree.X,inTree.Y,inTree.Z];
len = [0;len_tree(inTree)];
Idpar = ipar_tree(inTree);
pathLength = zeros(size(pointCloud,1),1);
for i = 1:size(pointCloud,1)
    [k,~] = dsearchn(B,pointCloud(i,:));
    pathLength(i) = sum(len(Idpar(k,:)+1));
end
end
