AlxSomaNeurites =[];
TransSomaNeurites = [];
DbxSomaNeurites = [];
BarhlSomaNeurites= [];

for i = 1:numel(cellIDs)
    if ismember(cellIDs{i},cellIDsAlx) ==1
        AlxSomaNeurites = [AlxSomaNeurites,length(allTrees{i}{1,1}{2})];
    elseif ismember(cellIDs{i},cellIDsTrans) ==1
        TransSomaNeurites = [TransSomaNeurites,length(allTrees{i}{1,1}{2})];
    elseif ismember(cellIDs{i},cellIDsDbx) ==1
        DbxSomaNeurites = [DbxSomaNeurites,length(allTrees{i}{1,1}{2})];
    else
        BarhlSomaNeurites = [BarhlSomaNeurites,length(allTrees{i}{1,1}{2})];
    end
end

SomaNeurites = [mean(AlxSomaNeurites);mean(TransSomaNeurites);mean(DbxSomaNeurites);mean(BarhlSomaNeurites)];
SomaNeuriteError = [std(AlxSomaNeurites);std(TransSomaNeurites);std(DbxSomaNeurites);std(BarhlSomaNeurites)];
