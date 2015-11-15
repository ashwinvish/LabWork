function [ tree ] = treeVisAV( swcPath )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

[tree, rawLength,SpecialPositions]= generateIrreducibleDoubleLinkedTree_WithDim(swcPath,[-1,0,2,3,4],2,true,[0.36,0.36,2]);
validNodes = [1:numel(tree)]; colorString = 'red'; newFigure = true; inducingNodes = [2]; highlightedNodes = [1];
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
                plot3(tempy,tempx,tempz,'black');
            end
        end
    end
end

for kk=1:numel(highlightedNodes)
    plot3(tree{highlightedNodes(kk)}{3}(2),tree{highlightedNodes(kk)}{3}(1),-tree{highlightedNodes(kk)}{3}(3),'Marker','o','MarkerSize',5,'Color','r')
end
for kk=1:size(SpecialPositions,1)
    plot3(SpecialPositions(kk,2),SpecialPositions(kk,1),-SpecialPositions(kk,3),'blue', 'LineWidth',2);
end

% for ll=1:numel(tree)
%     [x,y,z] =  [tree{1,ll}{3}(1),tree{1,ll}{3}(2),tree{1,ll}{3}(3)];
%     scatter3(x,y,-z,'Marker','o','MarkerSize','5','MarkerFaceColor','b');
% end


daspect([1,1,1]);
end

