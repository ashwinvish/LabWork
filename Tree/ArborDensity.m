val = nan(60,5);
bin = nan(60,5);
index =1;
for i = 1:numel(cellIDs) 
    if ismember (cellIDs{i}, cellIDsAlx)==1
        [d,x,y] = DendriticTree(allTrees{i},i,cellIDs,calx,true);
        val(1:length(x)+1,index) = [0,x]';
        bin(1:length(y),index) = y';
        index = index+1;
        clear x;
        clear y;
    else
        continue;
    end
end

figure;

for i = 1:numel(cellIDs)
    if ismember (cellIDs{i}, cellIDsTrans)==1
        DendriticTree(allTrees{i},i,cellIDs,ctrans,true);
    else
        continue;
    end
end
figure();
for i = 1:numel(cellIDs)
    if ismember (cellIDs{i}, cellIDsDbx)==1
        DendriticTree(allTrees{i},i,cellIDs,cdbx,true);
    else
        continue;
    end
end

figure();
for i = 1:numel(cellIDs)
    if ismember (cellIDs{i}, cellIDsL)==1
        DendriticTree(allTrees{i},i,cellIDs,cbarhl,true);
    else
        continue;
    end
end