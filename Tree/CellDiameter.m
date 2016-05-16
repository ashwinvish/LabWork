function [Dia, plengthAxon, plengthDen, AxnDia, DenDia]= CellDiameter( cellNo, allTrees, cellIDs, Display )
%CELLDIAMATER plots the cell of interest overlaid with the diameter at
%random locations
%   cellNo is the ID of the cell
%   allTrees is the cell with all the trees
%   cellIDs is the cell with the IDs of all the trees

% load stuff
load CellAxons.mat;
load CellDiameters.mat;


Diameter = Diameter{cellNo};
Dia(:,1) = 5*Diameter(:,1);
Dia(:,2) = 5*Diameter(:,2);
Dia(:,3) = -45*Diameter(:,3);

% display diameters on the cells
if Display == true
    figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(1,2,1);
    DisplayTree(allTrees{cellNo},[1],false, [eval([cellIDs{cellNo},'_axon'])],[0.8,0.8,0.8]);
    hold on;
    scatter3(5*Diameter(:,1), 5*Diameter(:,2), -45*Diameter(:,3), 50, Diameter(:,4),'filled');
    colormap jet;
    c = colorbar;
    c.FontSize = 40;
    c.Location = 'southOutside'; 
    c.Label.String = 'Diameter (nm)';
end


% get pathlenths to dendritic nodes
fname = sprintf('%s_WithTags.swc',cellIDs{cellNo});
[plength] = findPathLength(fname,allTrees{cellNo},[5,5,45],[5*Diameter(:,1),5*Diameter(:,2),45*Diameter(:,3)]);
%plength =  findPathLength_old(fname,[5,5,45],[5*Diameter(:,1),5*Diameter(:,2),45*Diameter(:,3)]);

% get identity of diameter node, axon or dendrite and the respective
% diameter

axnTree = AxonalTree(allTrees{cellNo},cellNo, cellIDs, false);
AxnDia = [];
DenDia = [];
index1 =1;
index2 = 1;
for i = 1:size(Dia,1)
    if ismember(Dia(i,:), axnTree) ==1
        AxnDia = [AxnDia;Diameter(i,:)]; % axonal node and diameter
        plengthAxon(index1) = plength(i);
        index1 = index1+1;
    else
        DenDia = [DenDia;Diameter(i,:)]; % dendritic node and diameter
        plengthDen(index2) = plength(i);
        index2 = index2+1;
    end
end
% normalize pathlenght


if isempty(axnTree)
    plengthAxon = 0;
    AxnDia = zeros(1,4);
end


plengthDen = plengthDen/max(plengthDen);
plengthAxon = plengthAxon/max(plengthAxon);



if Display == true
    subplot(1,2,2); %plot dendritic diameter and pathlength
    plot(plengthAxon/1000, AxnDia(:,4)./1000,'o', 'MarkerFaceColor','none', 'MarkerEdgeColor',[0,0.8,0], 'MarkerSize', 15, 'LineWidth',2);
    hold on;
    %plot(plengthAxon/1000, sort(axnFit/1000),'-','Color',[0,0.8,0]);
    plot(plengthDen/1000, DenDia(:,4)./1000,'o', 'MarkerFaceColor','none', 'MarkerEdgeColor',[0.9,0,0], 'MarkerSize', 15, 'LineWidth',2);
    xlabel('Norm. pathlength (\mum)', 'FontName', 'Arial', 'FontSize', 40);
    ylabel('Diameter (\mum)', 'FontName', 'Arial', 'FontSize', 40);
    set(gca,'LineWidth',2,  'FontName', 'Arial', 'FontSize', 40);
    legend({'Axon','Dendrite'});
    axis square;
    box off
end




end



