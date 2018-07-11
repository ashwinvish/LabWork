function [p1,p2,p3] = PlaneFit(PointCloud, PointColor, PlaneColor, ~)
% PLANEFIT is used to plot a plane over a cloud of points using the
% AFFINE_FIT function.
% PointCloud is an Nx3 matrix with the cartesian coordiated of the points
% PointColor is 3x1 vector of the rgb colors of the points.
% PlaneColor is a 3x1 vector of the rgb colors of the plane.
% ~ if nargin>3, Display a normal vector to the plane.


[p1,p2,p3] = affine_fit([PointCloud(:,1),PointCloud(:,2),-1*PointCloud(:,3)]);

XMaxMin = [max(PointCloud(:,1)),min(PointCloud(:,1))];
YMaxMin = [max(PointCloud(:,2)),min(PointCloud(:,2))];

Xgrid = linspace(XMaxMin(1), XMaxMin(2),5);
Ygrid= linspace(YMaxMin(1), YMaxMin(2),5);

[X,Y] = meshgrid(Xgrid,Ygrid);

% plot surface
surf(X,Y, -(p1(1)/p1(3)*X + p1(2)/p1(3)*Y-dot(p1,p3)/p1(3)),'FaceColor',PlaneColor,'FaceAlpha',0.5, 'FaceLighting','gouraud','EdgeColor', 'k', 'LineWidth', 2 );
hold on;

% plots point cloud

plot3(PointCloud(:,1),PointCloud(:,2), -PointCloud(:,3),'o','MarkerSize',10,'MarkerFaceColor',PointColor, 'MarkerEdgeColor','w');



    quiver3(p3(1),p3(2),p3(3), p1(1)/3, p1(2)/3,p1(3)/3, 50000,'k', 'LineWidth', 4);


%set orienataiton of plane and points
box on;
axis([ 20000 140000 60000 250000 -60000 0]);
daspect([1 1 1]); % make aspect ratio [1 1 1]
set (gca,'Ydir','reverse');
set (gca,'XTick',[], 'YTick',[],'ZTick', [], 'Ydir','reverse');
view([-180,90]);


end



