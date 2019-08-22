function [ logic ] =  isIBN(cellID)

% Check is the cell in cellID is an IBN cell
% logic = 1 if IBN cell
% logic = 0 if not an IBN cell
% type is the type of IBN cell

IBNall = [ 77125 77128 77135 77153 77231 77247 78685 77941 77137 79053 77940 ...
    77942 80971 77157 78550 79084 78557 79083 78567];

[logic] = ismember(cellID, IBNall);

end