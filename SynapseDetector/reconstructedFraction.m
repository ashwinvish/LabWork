function [fraction] = reconstructedFraction(cellID,df)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


for i = 1:numel(cellID)
    temp = SynapticPartners(cellID(i),1,df);
    A(i) = length(temp);
    temp = temp(temp<1e5);
    B(i) = length(temp);
    clear temp;
end

fraction = B./A;
end

