function [temp,AngleInDegrees,distanceToPlaneInMicrons] =  isMidline(cellID,isplot)

Midline

fname = '/usr/people/ashwinv/seungmount/research/Ashwin/Z-Brain/LowEMtoHighEM/SWC_all/consensus-20180412/';

if ischar(cellID)
    filename = sprintf('%s.swc',cellID);
else
    filename = sprintf('%d.swc',cellID);
end

swc = dlmread(fullfile(fname,filename), ' ');
coord = swc(:, 3:5) ;

if size(coord,1) < 250
    temp = NaN;
    AngleInDegrees = NaN;
    distanceToPlaneInMicrons = NaN;
    return;
else
    midPoints = [midlineDorsal; midlineVentral];
    midPoints(:,2) = midPoints(:,2);
    midPoints = midPoints.* [5,5,45];
    [p1,p2,p3] = affine_fit(midPoints);
    
    
    XMaxMin = [max(midPoints(:,1)),min(midPoints(:,1))];
    YMaxMin = [max(midPoints(:,2)),min(midPoints(:,2))];
    
    Xgrid = linspace(XMaxMin(1), XMaxMin(2),5);
    Ygrid= linspace(YMaxMin(1), YMaxMin(2),5);
    
    [X,Y] = meshgrid(Xgrid,Ygrid);
    
    vector1 = p1;
    % calculate vector to determine direction of the axon
    [topcoords,topIndex] = sort(coord(:,2),'ascend');
    %topcoords = sortedY(1:250);
    %topIndex = I(1:250);
    firstNonZerosIndex = [coord(topIndex(2:end),1)-coord(topIndex(1),1)~=0, ...
        coord(topIndex(2:end),2)-coord(topIndex(1),2)~=0, ...
        coord(topIndex(2:end),3)-coord(topIndex(1),3)~=0];
    
    firstNonZero = find((firstNonZerosIndex(:,1)+firstNonZerosIndex(:,2)+firstNonZerosIndex(:,3))==3);
    topNonZero = topIndex(firstNonZero(1)+1);
    
    secondNonZerosIndex = [coord(topIndex(firstNonZero(1)+2:end),1)-coord(topIndex(firstNonZero(1)+1),1)~=0, ...
        coord(topIndex(firstNonZero(1)+2:end),2)-coord(topIndex(firstNonZero(1)+1),2)~=0, ...
        coord(topIndex(firstNonZero(1)+2:end),3)-coord(topIndex(firstNonZero(1)+1),3)~=0];
    
    secondsNonZero =   find((secondNonZerosIndex(:,1)+secondNonZerosIndex(:,2)+secondNonZerosIndex(:,3))==3);
    
    %top1 = find(coord(:,2) == topcoords(1));
    top1 = topIndex(firstNonZero(1)+1+secondsNonZero(1));
    %top1 = topIndex(firstNonZero(1)+2);
    
    
    vector2 = coord(topNonZero,:)-coord(top1(1),:);
    vector2 = vector2 ./ norm(vector2);
    
    AngleInDegrees = atan2d(norm(cross(vector1,vector2)),dot(vector1,vector2));
    distanceToPlaneInMicrons = pdist2(p3,coord(topIndex(1),:))/1000;
    
    if (AngleInDegrees<40 && distanceToPlaneInMicrons < 50)
        temp = 1;
    else
        temp = 0;
    end
    
    % crossedPoints = find(coord(:,2)<p3(2));
    %
    % if crossedPoints>0
    %     temp = 1;
    % else
    %     temp = 0;
    % end
    
    
    if isplot
        
        scatter3(coord(:,1),coord(:,2),coord(:,3),'.');
        hold on;
        
        
        % plot surface
        surf(X,Y, -(p1(1)/p1(3)*X + p1(2)/p1(3)*Y-dot(p1,p3)/p1(3)),'FaceColor',[0.5,0.5,0.5],'FaceAlpha',0.5, 'FaceLighting','gouraud','EdgeColor', 'k', 'LineWidth', 2 );
        hold on;
        quiver3(p3(1),p3(2),p3(3), p1(1)/3, p1(2)/3,p1(3)/3, 50000,'B', 'LineWidth', 4);
        quiver3(coord(topNonZero,1),coord(topNonZero,2),coord(topNonZero,3), vector2(1),vector2(2),vector2(3), 50000,'R', 'LineWidth', 4);
        plot3(coord(1,1),coord(1,2),coord(1,3),'ko');
        
        set(gca,'XDir','reverse');
        
        
        xlabel('X');
        ylabel('Y');
        zlabel('Z');
        
    end
end
end

