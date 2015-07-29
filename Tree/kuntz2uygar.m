function [ tree ] = kuntz2uygar( filename,resolution)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

tree = load_tree(filename);
% tempx = tree.Y;
% tempy = tree.X;
% 
% tree.X = tempy;
% tree.Y = tempx;

tree.X = resolution(1)*tree.X;
tree.Y = resolution(2)*tree.Y;
tree.Z = -resolution(3)*tree.Z;

end

