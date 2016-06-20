function [TreeNodePoints] = NodePoints(tree)
x = zeros(numel(tree),1);
y = zeros(numel(tree),1);
z = zeros(numel(tree),1);
hold on;
for ll=1:numel(tree)
    tempx = tree{1,ll}{3}(1);
    tempy = tree{1,ll}{3}(2);
    tempz = tree{1,ll}{3}(3);
    x(ll,1) = tempx;
    y(ll,1) = tempy;
    z(ll,1) = tempz;
    %scatter3(tempy,tempx,-tempz,'Marker','o','MarkerFaceColor','b');
end
TreeNodePoints = [x,y,z];
end