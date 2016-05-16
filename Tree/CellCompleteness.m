figure();
plot(cellfun(@length, allPreSynapse),'ko' );
text(1:22,cellfun(@length, allPreSynapse),cellIDs);

figure();


%Display completed cells
figure();
index =1;
for i = 1:numel(cellIDs)
    if ismember(cellIDs{i}, cellIDsAlx) ==1
        subplot(4,2,index);
        %DisplayTree(allTrees{i},1,false,eval([cellIDs{i},'_axon']), [rand, rand, rand],allPreSynapse{i}, allPost{i});
        DisplayTree(allTrees{i},1,false,eval([cellIDs{i},'_axon']), [rand, rand, rand]);
        scatter3(allEnds{1,i}(:,1), allEnds{1,i}(:,2), -1*allEnds{1,i}(:,3), 'ko')
        title(cellIDs{i});
        view(-90,0);
        index = index +1;
    elseif  ismember(cellIDs{i}, cellIDsTrans) ==1
        subplot(4,2,index);
        %DisplayTree(allTrees{i},1,false,eval([cellIDs{i},'_axon']), [rand, rand, rand],allPreSynapse{i}, allPost{i});
         DisplayTree(allTrees{i},1,false,eval([cellIDs{i},'_axon']), [rand, rand, rand]);
         scatter3(allEnds{1,i}(:,1), allEnds{1,i}(:,2), -1*allEnds{1,i}(:,3), 'ko')
        title(cellIDs{i});
        view(-90,0);
        index = index +1;
        
    end
end

% complete cells are Int1_4, Int1_5 and Int1_6

completeLength = [allRawLength{4},allRawLength{5},allRawLength{6}];
completePreSynapses = [length(allPreSynapse{4}), length(allPreSynapse{5}), length(allPreSynapse{6})];
completePostSynapses = [length(allPost{4}), length(allPost{5}), length(allPost{6})];
