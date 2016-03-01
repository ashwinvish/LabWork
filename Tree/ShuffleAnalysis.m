% Shuffle analysis
clc;

ShuffleSteps = 1000*[1,5,10,15,20];

for i = 1:numel(ShuffleSteps)
    PotSites  = Shuffle(allTrees,cellIDs,cellIDsAxon,false,ShuffleSteps(i),1000);
    ShuffName = sprintf('Shuffle%dum_%s.mat',ShuffleSteps(i),date);
    save(ShuffName,'PotSites');
end
