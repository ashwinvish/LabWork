function [dist,pair]=PairwiseSpatialStructure2(R,M,per)
%R is a Nx3 matrix of cellular positions, the first column is the rc
%position, the second column is the ML position and the third column is the
%DV position, i.e. R=[X Y Z];
%RC spatial structure
%if per is 1 you compute pairwise percent difference otherwise you compute
%absolute difference
% R(:,2)=abs(R(:,2));
for q=1:3
    X=R(:,q);
    if q==2;X=abs(X);end
    DX=[];DM=[];
    for k=1:2;
        if k==1;a=find(R(:,2)>0);else;a=find(R(:,2)<0);end
        yy=X(a);mns=M(a);
        clear dx dm
        if length(a)>0
            for i=1:length(a);
                for j=1:length(a);
                    dx(i,j)=yy(i)-yy(j);
                    if per==1;
                        dm(i,j)=(mns(i)-mns(j))/(mns(i)+mns(j))/2;
                    else
                        dm(i,j)=mns(i)-mns(j);
                    end
                end
            end
            a=find(eye(size(dx))==1);dx(a)=[];dm(a)=[];
            DX=[DX dx];DM=[DM dm];
        end
    end
    %     [a,b]=sort(DX);DX=DX(b)';DM=DM(b)';
    dist(:,q)=DX(:);pair(:,q)=DM(:);
end