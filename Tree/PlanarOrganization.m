% Organization of AlxPost and DbxPost sites

P1 = PlaneFit(AlxPost,calx,calx);
hold on;
P2 = PlaneFit(DbxPost, cdbx, cdbx);
view(-90,0);

AlxDbxPostAngle = acosd(dot(P1,P2));
if AlxDbxPostAngle>90
    AlxDbxPostAngle = 180-AlxDbxPostAngle;
end

str1 = sprintf('Angle between AlxPost and DbxPost plane is %3.2fd', AlxDbxPostAngle);
title(str1);

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
str2 = sprintf('Angle between AlxPost and DbxPost plane is %3.2fd', AlxPrePostAngle);
title(str2);

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

str3 = sprintf('Angle between AlxPost and DbxPost plane is %3.2f', AlxPreDbxPostAngle);
title(str3);

% All postSynaptic planes
figure();

PlaneFit(AlxPost,calx,calx);
hold on;
PlaneFit(DbxPost, cdbx, cdbx);
PlaneFit(TransPost, ctrans, ctrans);
PlaneFit(BarhlPost, cbarhl, cbarhl);

view(-90,0);
