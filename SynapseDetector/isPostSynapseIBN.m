function [logic] = isPostSynapseIBN(cellID,df)
%UNTITLED14 Summary of this function goes here
%   Detailed explanation goes here

IBNall = [82714,82710,79083,82264,78550,79084,82712,82709,82713,83184,...
    77125,77942,77231,77247,77153,78685,77137,83183];

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

