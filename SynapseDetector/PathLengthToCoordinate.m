function [pathLength] = PathLengthToCoordinate(A,inTree)
B = [inTree.X,inTree.Y,inTree.Z];
len = [0;len_tree(inTree)];
Idpar = ipar_tree(inTree);

pathLength = zeros(size(A,1),1);
for i = 1:size(A,1)
    [k,~] = dsearchn(B,A(i,:));
    pathLength(i) = sum(len(Idpar(k,:)+1));
end
end
