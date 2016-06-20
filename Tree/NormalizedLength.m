function [NormPathLength] = NormalizedLength( PathLength )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

NormPathLength = PathLength/max(PathLength);

end

