% plot Trees

col = distinguishable_colors(24);
DisplayTree(allTrees{5},[],true,eval([cellIDs{5},'_axon']),col(5,:))                           % 1_5
hold on;
TreeSomata(5,col(5,:));
DisplayTree(allTrees{2},[],false,eval([cellIDs{2},'_axon']),col(2,:));                         % 1_2
TreeSomata(2,col(2,:));
DisplayTree(allTrees{9},[],false,eval([cellIDs{9},'_axon']),col(9,:))    ;                     % 2_2
TreeSomata(9,col(9,:));
DisplayTree(allTrees{22},[],false,eval([cellIDs{22},'_axon']), col(22,:))                      % 3_6
TreeSomata(22,col(22,:));
DisplayTree(allTrees{12},[],false,eval([cellIDs{12},'_axon']), col(12,:))                      % 2_5
TreeSomata(12,col(12,:));

hold on;

% plot actual synapses

syn1 = [16662	24415	470
    18872	23363	468
    19645	36619	960
    19279	35061	1046];
syn1(:,3) = -1*syn1(:,3);
s1 = scatter3(5*syn1(:,1),5*syn1(:,2),45*syn1(:,3),'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','none');
s1Markers = s1.MarkerHandle;
s1.SizeData = 100;


syn2 = [14363 21117	390
    16022 21879	364];
syn2(:,3) = -1*syn2(:,3);

s2 = scatter3(5*syn2(:,1),5*syn2(:,2),45*syn2(:,3),'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','none');
s2Markers = s2.MarkerHandle;
s2.SizeData = 100;

h1 = gcf;
set(h1,'color','none', 'units','normalized','outerposition',[0 0 1 1]);
axis vis3d;
PlotViews(h1);

 
 