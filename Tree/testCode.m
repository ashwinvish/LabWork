figure();
index =1;
for i = 1:length(cellIDs)
    if ismember(cellIDs(i),cellIDsAlx) == 1
        [A,B,C,D,E] = CellDiameter(i,allTrees, cellIDs, false);
        if isempty(D)
            AxnDia{i} = [];
            DenDia{i} = E(:,4);
        else
            AxnDia{i} = D(:,4);
            DenDia{i} = E(:,4);
            AxnPlength = B;
            DenPlength = C;
        end
        AlxAxnDia(index,1:size(AxnDia{i},1)) = AxnDia{i}';
        AlxDenDia(index,1:size(DenDia{i},1)) = DenDia{i}';
        AlxAxonPlength(index,1:size(AxnPlength,1)) = AxnPlength;
      %  AlxDenPlength(index,1:size(DenPlength,1)) = DenPlength;
        plot([1,3],[mean(AxnDia{i})/1000,mean(DenDia{i})/1000], '-o','MarkerFaceColor',[0.8,0.8,0.8],'MarkerSize', 35,'MarkerEdgeColor','k','Color','k', 'LineWidth',2);
        hold all;
        index = index+1;
    end
end