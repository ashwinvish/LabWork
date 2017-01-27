AlxSomaTime = [];
TransSomaTime = [];
DbxSomaTime  = [];
BarhlSomaTime = [];
CellColor = [];

for i = 1:numel(cellIDs)
    if ismember(cellIDs{i},cellIDsAlx) ==1
        AlxSomaTime = [AlxSomaTime; tau(i)];
    elseif ismember(cellIDs{i}, cellIDsTrans) ==1
        TransSomaTime = [TransSomaTime; tau(i)];
    elseif ismember(cellIDs{i}, cellIDsDbx) ==1
        DbxSomaTime = [DbxSomaTime; tau(i)];
    else
        BarhlSomaTime = [BarhlSomaTime; tau(i)];
    end
end
