function [dist,rhodist] = PlotPairs( cell1, cell2 , rho1, rho2, rho, map, varargin)
%   PlotPairs function to plot a pair of cells
%   cell1 is the ID of the first cell
%   cell2 is the ID of the second cell
%   rho1 is the presistence measure for cell1
%   rho2 is the persistence measure for cell2
%   rho is the persistence measuer for all cells
%   varargin ture or false to plot synapses

cells = {cell1, cell2};
cmap1 = find(rho1 == sort(rho),1);
cmap2 = find(rho2 == sort(rho),1);

%cmap1 = find(rho1 == rho);
%cmap2 = find(rho2 == rho);

if nargin >6 % if you want to plot synapse locations
    
    temp1 = [allPreSynapse{cell1},allPreSynapse{cell2}] ;
    temp2 = [allPostSynapse{cell1},allPostSynapse{cell2}];
    
    treeVisualizer(cells{1}, [1],[],[temp1(1) temp2(1)],false,{map(cmap1,:)}, 1:numel(cells{1}), false);
    treeVisualizer(cells{2}, [1],[],[temp1(2) temp2(2)],false,{map(cmap2,:)}, 1:numel(cells{2}), false);
else
    treeVisualizer(cells{1}, [1],[],[],false,{map(cmap1,:)}, 1:numel(cells{1}), false);
    treeVisualizer(cells{2}, [1],[],[],false,{map(cmap2,:)}, 1:numel(cells{2}), false);
end

cell1Soma = cell1{1}{1,3};
cell2Soma = cell2{1}{1,3};

dist = sqrt((cell1Soma(1)-cell2Soma(1))^2 +(cell1Soma(2)-cell2Soma(2))^2 +(cell1Soma(3)-cell2Soma(3))^2 );
rhodist = abs(rho1-rho2);

%title(sprintf('%s and %s \n rho1: %.1f and rho2: %.1f', cellIDs{I(i,1)}, cellIDs{I(i,2)}, rho(I(i,1)), rho(I(i,2)) ));

end

