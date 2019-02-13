function [value] = isContra( cellId )
% isContra query to determine is cellId is a contra cell

if ismac
    contraList = dlmread('/Users/ashwin/Documents/SynapseDetector/contraList-01242019.txt','\t',1,0);
else
    contraList = dlmread('/usr/people/ashwinv/seungmount/research/Ashwin/SynapseDetector/contraList-07032018.txt','\t',1,0);
end

%location = ismember(cellId,contraList);

for i = 1:size(cellId,1)
         l = find(cellId(i) == contraList(:,1));
    if ~isempty(l)
        if isnumeric(contraList(l,3))
            value(i) = logical(contraList(l,3));
        else
            value(i) = false;
        end
    else
        value(i) = false;
    end
end
    value = value';

end

