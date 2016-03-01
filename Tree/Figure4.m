col = distinguishable_colors(6,calx);
DisplayTree(allTrees{5},[1],false, [eval([cellIDs{5},'_axon'])],col(2,:),allPreSynapse{5}, allPostSynapse{5});
DisplayTree(allTrees{6},[1],false, [eval([cellIDs{6},'_axon'])],col(3,:),allPreSynapse{6}, allPostSynapse{6});

axis vis3d;
set(gca, 'Color', [calx,0.2]);
set(gcf,'Color','none');

%%

col = distinguishable_colors(2,ctrans);
DisplayTree(allTrees{7},[1],false, [eval([cellIDs{7},'_axon'])],col(1,:),allPreSynapse{7}, allPostSynapse{7});
DisplayTree(allTrees{21},[1],false, [eval([cellIDs{21},'_axon'])],col(2,:),allPreSynapse{21}, allPostSynapse{21});

axis vis3d;
set(gca, 'Color', [ctrans,0.2]);
set(gcf,'Color','none');

%%

col = distinguishable_colors(8,cdbx);
DisplayTree(allTrees{2},[1],false, [eval([cellIDs{2},'_axon'])],col(1,:),allPreSynapse{2}, allPostSynapse{2});
DisplayTree(allTrees{9},[1],false, [eval([cellIDs{9},'_axon'])],col(4,:),allPreSynapse{9}, allPostSynapse{9});

axis vis3d;
set(gca, 'Color', [cdbx,0.2]);
set(gcf,'Color','none');

%%

col = distinguishable_colors(6,cbarhl);
DisplayTree(allTrees{1},[1],false, [eval([cellIDs{1},'_axon'])],col(1,:),allPreSynapse{1}, allPostSynapse{1});
DisplayTree(allTrees{20},[1],false, [eval([cellIDs{20},'_axon'])],col(6,:),allPreSynapse{20}, allPostSynapse{20});

axis vis3d;
set(gca, 'Color', [cbarhl,0.2]);
set(gcf,'Color','none');