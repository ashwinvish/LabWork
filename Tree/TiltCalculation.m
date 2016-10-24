
MauthnerCell = [5*14474,5*49530,-45*448];                                                       % cartesian coordinates for the center of the Mauthner cell
Mid2C1 = [5*9344,5*44812,-45*468];
Mid2C2 = [5*10816,5*45268,-45*466];
Mid3C1 = [5*13944,5*42036,-45*137];
Mid3C2 = [5*13072,5*40892,-45*28];
CAD = [5*6846,5*36655, -45*1 ];
CAV = [5*6894,5*38023,-45*16];

MeanMid2 = mean([Mid2C1;Mid2C2]);
MeanMid3 = mean([Mid3C1;Mid3C2]);
MeanCa = mean([CAD;CAV]);

AllMeans = [MauthnerCell;MeanMid2;MeanMid3;MeanCa];
[p1,p2,p3] = PlaneFit([AllMeans(:,1), AllMeans(:,2), -1*AllMeans(:,3)], [1,0,0],[0,1,0],1);
close all;
%Find midline between Ca and Mid3

DisplayTree(allTrees{6},[],false, [eval([cellIDs{6},'_axon'])],[0.1,0.8,0.1], allPreSynapse{6}, allPost{6});
TreeSomata(6,[0.1,0.8,0.1]);
DisplayTree(allTrees{4},[],false, [eval([cellIDs{4},'_axon'])],[0.8,0.1,0.1], allPreSynapse{4}, allPost{4});
TreeSomata(4,[0.8,0.1,0.1]);

Border6_7 = mean([MeanMid3; MeanCa]);
Border6_7_projecton = Projection(Border6_7, p1,p3);


index = 2;
Z = 0:1000:60000;
X = Border6_7_projecton(1)* ones(1,length(Z));
Y(1) = Border6_7_projecton(2);
for i = 2:length(Z)
    Y(index) = Y(index-1)-370; % 370/292 nm
    index = index+1;
end
hold on;
scatter3(X,Y,-1*Z,300, 'k.');
scatter3(MauthnerCell(1), MauthnerCell(2), MauthnerCell(3),500,'k', 'p', 'MarkerFaceColor','k');
scatter3(MeanMid2(1), MeanMid2(2), MeanMid2(3),300,'r', 'p', 'MarkerFaceColor','r', 'MarkerEdgeColor','k');
scatter3(MeanMid3(1), MeanMid3(2), MeanMid3(3),300,'b','p', 'MarkerFaceColor','b', 'MarkerEdgeColor','k');
scatter3(MeanCa(1), MeanCa(2), MeanCa(3),300, 'g','p', 'MarkerFaceColor','g', 'MarkerEdgeColor','k');


