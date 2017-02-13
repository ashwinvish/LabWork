% CellGroupColor.m record the colors for the cells by group 
CellColor = [];
for i = 1:length(cellIDs)
    if ismember(cellIDs(i), cellIDsAlx) ==1
        CellColor(i,:) = calx;
    elseif ismember(cellIDs(i), cellIDsTrans) ==1
        CellColor(i,:) = ctrans;
    elseif ismember(cellIDs(i),cellIDsDbx) == 1
        CellColor(i,:) = cdbx;
    else
        CellColor(i,:) = cbarhl;
    end
end
