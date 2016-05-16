function [sh,h] = TreeSomata( treeNo, col )
%TREESOMATA plots the 3d shape for the somata of the tree
%   treeNo is the serial number of the tree, e.g. In1_1 is 1
%   col is the color of the shape object

load CellSomata.mat;
res = [5,5,-45]; % convert from pixels to nm

sh = alphaShape(res(1)*TreeSoma{treeNo}(:,1), res(2)*TreeSoma{treeNo}(:,2), res(3)*TreeSoma{treeNo}(:,3));
sh.Alpha;
sh.Alpha = sh.Alpha + 1000;
h = plot(sh);
lightangle(-155,30);
h.FaceColor =  col;
h.EdgeColor = 'none';
h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.6;
h.DiffuseStrength = 0.1;
h.SpecularStrength = 0.2;
h.SpecularExponent = 25;
h.BackFaceLighting = 'reverselit';

axis([ 20000 140000 60000 250000 -60000 0]);
daspect([1 1 1]);
end

