
index = 1;
for i = 1:numel(cellIDs)
    if ismember(cellIDs{i},cellIDsDbx)==1
            DbxPostSynapticPoints{index} = allPost{i};
            DbxSomata(index,:) = CellSoma(i,:);
            index = index+1;
    end
end
A = [];
translatedSomata(1,:) = DbxSomata(1,:);
translatedPostSynapticPoints{1} = DbxPostSynapticPoints{1};
for i = 2:length(DbxPostSynapticPoints)
    [translatedPostSynapticPoints{i},translatedSomata(i,:)] = translateTree(DbxPostSynapticPoints{i},DbxSomata(i,:), translatedSomata(1,:));
    A = [A ;translatedPostSynapticPoints{i}];
end

col =  cbrewer('qual','Set1' , 9);

    for i = 1:numel(translatedPostSynapticPoints)
        plot3(translatedSomata(i,1),translatedSomata(i,2), translatedSomata(i,3),'o','color',col(i,:))
        scatter3(translatedPostSynapticPoints{i}(:,1), translatedPostSynapticPoints{i}(:,2), translatedPostSynapticPoints{i}(:,3),50,...
            'MarkerEdgeColor', col(i,:), 'MarkerFaceColor', col(i,:));
        hold on;
    end
    
    view(gca, -180,0);
    axis(gca, 'normal');
    set (gca,'XTick',[], 'YTick',[],'ZTick', [],'XColor','none','YColor','none', 'ZColor','none');
    box(gca,'off');

%%
figure();
P1 = PlaneFit(A,col(1,:),col(1,:));
hold on;
P2 = PlaneFit(DbxPost, cdbx, cdbx);
view(-90,0);

AlxDbxPostAngle = acosd(dot(P1,P2));
if AlxDbxPostAngle>90
    AlxDbxPostAngle = 180-AlxDbxPostAngle;
end

% str1 = sprintf('Angle between AlxPost and DbxPost plane is %3.2fd', AlxDbxPostAngle);
% title(str1);
% hold on;
% scatter3(MauthnerCell(1,1),MauthnerCell(1,2),MauthnerCell(1,3), 1000,'p','MarkerFaceColor','k', 'MarkerEdgeColor', 'k');
% scatter3(Mid2C1(1,1),Mid2C1(1,2),Mid2C1(1,3), 500,'p','MarkerFaceColor','r', 'MarkerEdgeColor', 'k');
% scatter3(Mid2C2(1,1),Mid2C2(1,2),Mid2C2(1,3), 500,'p','MarkerFaceColor','r', 'MarkerEdgeColor', 'k');
% scatter3(Mid3C1(1,1),Mid3C1(1,2),Mid3C1(1,3), 500,'p','MarkerFaceColor','b', 'MarkerEdgeColor', 'k');
% scatter3(Mid3C2(1,1),Mid3C2(1,2),Mid3C2(1,3), 500,'p','MarkerFaceColor','b', 'MarkerEdgeColor', 'k');
% scatter3(CAD(1,1), CAD(1,2), CAD(1,3), 500,'p','MarkerFaceColor','g', 'MarkerEdgeColor', 'k');
% scatter3(CAV(1,1), CAV(1,2), CAV(1,3), 500,'p','MarkerFaceColor','g', 'MarkerEdgeColor', 'k');


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