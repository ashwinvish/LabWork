
% plot Trees

col = distinguishable_colors(24);
DisplayTree(allTrees{5},[],true,eval([cellIDs{5},'_axon']),col(5,:))                           % 1_5
hold on;
TreeSomata(5,col(5,:));
DisplayTree(allTrees{2},[],false,eval([cellIDs{2},'_axon']),col(2,:));                     % 1_2
TreeSomata(2,col(2,:));
DisplayTree(allTrees{9},[],false,eval([cellIDs{9},'_axon']),col(9,:))    ;% 2_2
TreeSomata(9,col(9,:));
DisplayTree(allTrees{22},[],false,eval([cellIDs{22},'_axon']), col(22,:))                    % 3_6
TreeSomata(22,col(22,:));

hold on;

% plot actual synapses

 syn = [16662	24415	470
 18872	23363	468
 19645	36619	960
 19279	35061	1046];
 syn(:,3) = -1*syn(:,3);
 s2 = scatter3(5*syn(:,1),5*syn(:,2),45*syn(:,3),'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','none');
 s2Markers = s2.MarkerHandle;
 %s2Markers.FaceColorData = uint8(255*[1;0;0;0.3]);  % Alpha=0.3 => 70% transparent red
 s2.SizeData = 100;