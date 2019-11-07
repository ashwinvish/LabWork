function [logic] = isPostSynapseIBN(cellID,df)
%UNTITLED14 Summary of this function goes here
%   Detailed explanation goes here

IBNall = [ 77125 77128 77135 77153 77231 77247 78685 77941 77137 79053 ...
    77940 77942 80971 77157 78550 79084 78557 79083 78567];

logic = zeros(size(cellID));
for i = 1:numel(cellID)
    [A,~] = SynapticPartners(cellID(i),2,df);
    valid = ismember(A,IBNall);
    if sum(valid)>0
        logic(i) = 1;
    end
    clear A;
    clear valid;
end
end

