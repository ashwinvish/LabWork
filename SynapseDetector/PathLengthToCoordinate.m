function [pathLength] = PathLengthToCoordinate(A,inTree)

B = [inTree.X,inTree.Y,inTree.Z];
len = [0;len_tree(inTree)];
Idpar = ipar_tree(inTree);

pathLength = zeros(length(A),1);
for i = 1:length(A)
    [k,~] = dsearchn(B,A(i,:));
    pathLength(i) = sum(len(Idpar(k,:)+1));
end
end
