% Potentail Synapse Analysis

PotSites = [];
PotentialStites = [];
index = 1;

for i = 1:numel(cellIDs)
     if ismember(cellIDs(i),cellIDsAlx) == 1
         figure;
         for j = 1:numel(cellIDs)
         subplot(3,8,j); 
         [PotSitesTemp]= PotentialSites(allTrees, cellIDs, i,j,true,false);
         [PotSites] = [PotSitesTemp; PotSites];
         pause (1);
         end
         suptitle(cellIDs(i));
     end
     PotentialStites{index} = PotSites;
     index = index+1;
end