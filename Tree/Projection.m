function [proj] =  Projection(Point, unitVector,PointOnPlane)
% PROJECTION calculated the projection of a point on a plane
    proj = Point'-dot((Point-PointOnPlane),unitVector)*unitVector;
end
