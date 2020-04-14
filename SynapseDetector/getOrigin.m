
function [origin] = getOrigin(cellID)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

for i = 1:numel(cellID)
    if (isExistReRoot(cellID(i)) == 1)
        tree = SwctoZbrian(cellID(i));
        origin(i,:) = [tree{1}.X(1), tree{1}.Y(1), tree{1}.Z(1)];
        clear tree;
    else
        origin(i,:) = [0,0,0];
    end
end
end

