%SynapseLocations

AlxPostSites = []
AlxPreSites = [];
TransPostSites = [];
TransPreSites = [];
DbxPostSites = [];
BarhlPostSites = [];

for i = 1:numel(cellIDs)
    if ismember(cellIDs{i}, cellIDsAlx) ==1
        AlxPostSites = [AlxPostSites; allPost{i}];
        AlxPreSites = [AlxPreSites; allPreSynapse{i}];
    elseif ismember(cellIDs{i}, cellIDsTrans) ==1
        TransPostSites = [TransPostSites; allPost{i}];
        TransPreSites = [TransPreSites; allPreSynapse{i}];
    elseif ismember(cellIDs{i}, cellIDsDbx) ==1
        DbxPostSites = [DbxPostSites; allPost{i}];
    else 
        BarhlPostSites = [BarhlPostSites; allPost{i}];
    end
end

