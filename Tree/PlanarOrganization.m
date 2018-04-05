% Organization of AlxPost and DbxPost sites
figure();
P1 = PlaneFit(AlxPost,calx,calx);
hold on;
P2 = PlaneFit(DbxPost, cdbx, cdbx);
view(-90,0);

AlxDbxPostAngle = acosd(dot(P1,P2));
if AlxDbxPostAngle>90
    AlxDbxPostAngle = 180-AlxDbxPostAngle;
end

% str1 = sprintf('Angle between AlxPost and DbxPost plane is %3.2fd', AlxDbxPostAngle);
% title(str1);
hold on;
scatter3(MauthnerCell(1,1),MauthnerCell(1,2),MauthnerCell(1,3), 1000,'p','MarkerFaceColor',map(1,:), 'MarkerEdgeColor', 'k');
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


% Angle between AlxPre and AlxPost
figure();

P3 = PlaneFit(AlxPost,[0.8,0,0],[0.8,0,0]);
hold on;
P4 = PlaneFit(AlxPre,[0,0.9,0],[0,0.9,0]);
set(gca,'color', [calx,0.2]);
view (-90,0);
AlxPrePostAngle = acosd(dot(P3,P4));
if AlxPrePostAngle>90
    AlxPrePostAngle = 180-AlxPrePostAngle;
end
% str2 = sprintf('Angle between AlxPost and DbxPost plane is %3.2fd', AlxPrePostAngle);
% title(str2);

% Angle between AlxPre and DbxPost
figure();

P5 = PlaneFit(AlxPre,[0,0.9,0],calx);
hold on;
P6 = PlaneFit(DbxPost,[0.8,0,0],cdbx);
view(-90,0);
AlxPreDbxPostAngle = acosd(dot(P3,P4));
if AlxPreDbxPostAngle>90
    AlxPreDbxPostAngle = 180-AlxPreDbxPostAngle;
end

% str3 = sprintf('Angle between AlxPost and DbxPost plane is %3.2f', AlxPreDbxPostAngle);
% title(str3);

% All postSynaptic planes
figure();

PlaneFit(AlxPost,calx,calx);
hold on;
PlaneFit(DbxPost, cdbx, cdbx);
PlaneFit(TransPost, ctrans, ctrans);
PlaneFit(BarhlPost, cbarhl, cbarhl);

view(-90,0);
