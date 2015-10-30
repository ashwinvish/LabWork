function AnatomicalSpline( X , Y, Z , color)
%ANATOMICALSPLINES draws a preiodic interpolating spline through points
%   X is the X points in cartesian coordinates
%   Y is the Y points in cartesian coordinates
%   Z is the Z points in cartesian coordinates
%   COLOR [0.1 0.1 0.1]

zVals = unique(Z);
x = [];y = []; z = [];

for i = 1:numel(zVals)
    for j = 1:length(X)
        if zVals(i) == Z(j)
            [x] = [x ; X(j)];
            [y] = [y ; Y(j)];
            [z] = [z ; Z(j)];
        end
    end
    hold all;
    x = x';
    y = y';
    z = z';
    xyz = [x;y;z];
    fnplt(cscvn(xyz(:,[1:end 1])),'k',1);
    [x] = []; [y] = []; [z] = [];
end

end

