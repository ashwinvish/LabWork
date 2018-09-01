function [value] = isContra( cellId )
% isContra query to determine is cellId is a contra cell

contraList = dlmread('/usr/people/ashwinv/seungmount/research/Ashwin/SynapseDetector/contraList-07032018.txt','\t',1,0);
location = find(contraList(:,1) == cellId);
value = contraList(location,3);

end

