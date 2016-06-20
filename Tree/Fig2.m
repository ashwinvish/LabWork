
%1_6

DisplayTree(allTrees{6},[],false, [eval([cellIDs{6},'_axon'])],[0.3 0.3 0.3],allPreSynapse{6}, allPostSynapse{6});
hold on
TreeSomata(6,[0.3,0.3,0.3]);
box off
set(gca, 'color','none');
set(gca, 'XColor', 'none', 'YColor','none', 'ZColor', 'none');
BeginMyelin = [5*7493,5*40208,-45*725;5*7056,5*46560,-45*1200;5*7478, 5*37592,-45*531];
EndMyelin = [5*7366, 5*44816, -45*1111;5*7367, 5*48064, -45*1278;5*7871, 5*36291, -45*455];
scatter3(BeginMyelin(:,1), BeginMyelin(:,2), BeginMyelin(:,3), 'MarkerFaceColor', 'y')
scatter3(EndMyelin(:,1), EndMyelin(:,2), EndMyelin(:,3), 'MarkerFaceColor', 'y')

%1_5

DisplayTree(allTrees{5},[],false, [eval([cellIDs{5},'_axon'])],[0.3 0.3 0.3],allPreSynapse{5}, allPostSynapse{5});
hold on
TreeSomata(5,[0.3,0.3,0.3]);
set(gca, 'XColor', 'none', 'YColor','none', 'ZColor', 'none');
set(gca, 'color','none');

%2_2

DisplayTree(allTrees{9},[],false, [eval([cellIDs{9},'_axon'])],[0.3 0.3 0.3],allPreSynapse{9}, allPostSynapse{9});
hold on;
TreeSomata(9,[0.3,0.3,0.3]);
set(gca, 'color','none');
set(gca, 'XColor', 'none', 'YColor','none', 'ZColor', 'none');

%3_6
DisplayTree(allTrees{21},[],false, [eval([cellIDs{21},'_axon'])],[0.3 0.3 0.3],allPreSynapse{21}, allPostSynapse{21});
hold on
TreeSomata(21,[0.3,0.3,0.3]);
set(gca, 'XColor', 'none', 'YColor','none', 'ZColor', 'none');
set(gca, 'color','none');
