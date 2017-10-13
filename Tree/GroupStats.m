%GroupStats

AlxLength = [];
TransLength = [];
DbxLength = [];
BarhlLength = [];

AlxAxLength = [];
TransAxLength = [];
DbxAxLength = [];
BarhlAxLength = [];

AlxDenLength = [];
TransDenLength = [];
DbxDenLength = [];
BarhlDenLength = [];




for i = 1:numel(cellIDs)
    if ismember(cellIDs{i}, cellIDsAlx) ==1
        AlxLength = [AlxLength, allRawLength{i}];
        AlxAxLength = [AlxAxLength; axLength(i)];
        AlxDenLength = [AlxDenLength, denLength(i)];
        
    elseif ismember( cellIDs{i}, cellIDsTrans) ==1
        TransLength = [TransLength;allRawLength{i}];
        TransAxLength = [TransAxLength; axLength(i)];
        TransDenLength = [TransDenLength; denLength(i)];
        
    elseif ismember( cellIDs{i}, cellIDsDbx) ==1
        DbxLength = [DbxLength; allRawLength{i}];
        DbxAxLength = [DbxAxLength;axLength(i)];
        DbxDenLength = [DbxDenLength; denLength(i)];
    else
        BarhlLength = [BarhlLength; allRawLength{i}];
        BarhlAxLength = [BarhlAxLength;axLength(i)];
        BarhlDenLength = [BarhlDenLength; denLength(i)];
    end
end

