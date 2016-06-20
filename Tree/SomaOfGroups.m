% SomaByGroups
AlxSoma = [];
TransSoma = [];
DbxSoma = [];
BarhlSoma = [];

for i = 1:numel(cellIDs)
    if ismember(cellIDs{i},cellIDsAlx) ==1
        AlxSoma = [AlxSoma; CellSoma(i,:)];
    elseif ismember(cellIDs{i}, cellIDsTrans) ==1
        TransSoma = [TransSoma; CellSoma(i,:)];
    elseif ismember(cellIDs{i}, cellIDsDbx) ==1
        DbxSoma = [DbxSoma; CellSoma(i,:)];
    else
        BarhlSoma = [BarhlSoma; CellSoma(i,:)];
    end
end

