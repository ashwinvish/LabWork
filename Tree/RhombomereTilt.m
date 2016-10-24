% generatge figures for rhombomeres and Abducens

imInfo = imfinfo('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/Experiements/10122012-1/ZFishAllWafers_292nm_1748sec-cleanAlignment.tif');
imWidth = imInfo(1).Width;
imHeight = imInfo(1).Height;
imNumber = length(imInfo);
TifLink = Tiff('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/Experiements/10122012-1/ZFishAllWafers_292nm_1748sec-cleanAlignment.tif','r');
for i = 1:imNumber
    TifLink.setDirectory(i);
    FinalImage(:,:,i) = TifLink.read();
end

% anatomical landmarks
MCell = [210,	828,	458];
M2C1 = 	[117,	736,	533];
M2C2 = 	[150,	748,	398];
MeanM2C = mean([M2C1; M2C2]);
M3C1 = 	[189,	653,	47];
M3C2 = 	[196,	670,	115];
MeanM3C = mean([M3C1; M3C2]);
CaD =   [80,    583,    1];
CaV =   [86,    604,    23];
MeanCa  = mean([CaD; CaV]);

% Find plane through all points
figure();
PC = [MCell;MeanM2C; MeanM3C; MeanCa];
[p1,p2,p3 ]  = affine_fit([PC(:,1), PC(:,2), PC(:,3)]);

XMaxMin = [max(PC(:,1)),min(PC(:,1))];
YMaxMin = [max(PC(:,2)),min(PC(:,2))];

Xgrid = linspace(XMaxMin(1), XMaxMin(2),10);
Ygrid= linspace(YMaxMin(1), YMaxMin(2),10);

[X,Y] = meshgrid(Xgrid,Ygrid);

% plot surface
surf(X,Y, -(p1(1)/p1(3)*X + p1(2)/p1(3)*Y-dot(p1,p3)/p1(3)),'FaceColor','r','FaceAlpha',0.5, 'FaceLighting','gouraud','EdgeColor', 'k' );
hold on;
set(gca, 'XLim', [0, 511],'YLim', [0,903], 'ZLim', [0, 1747]);
% plots point cloud

plot3(PC(:,1),PC(:,2),PC(:,3),'o','MarkerSize',6,'MarkerFaceColor','r', 'MarkerEdgeColor','k');
%imshow(FinalImage(:,:,150));

clear X;
clear Y;

[X,Y] = meshgrid(linspace(0, 511, 511), linspace(0,903,903));           % image dimensions are 903x511
imageplane = [458, 1000]; % 1 to 1285

for i =  1: length(imageplane)
    Z = imageplane(i)*ones(size(X));
    hsurf = surface(X,Y,Z,histeq(FinalImage(:,:,imageplane(i)),120),'FaceColor','texturemap','EdgeColor','none');
    %alpha(hsurf, 1);
    colormap gray;
end
set(gca, 'YDir', 'reverse');
daspect([1,1,6.48]);

%Integrator Cells

Small(1,:) = [382	482	1];
Small(2,:) = [195	366	65];
Small(3,:) = [154	322	22];
Small(4,:) = [100	614	70];
Small(5,:) = [108	640	14];
Small(6,:) = [174	643	14];
Small(7,:) = [170	613	14];
Small(8,:) =  [208	396	250];
Small(9,:) = [192	350	161];
Small(10,:) = [180 353 250];
Small(11,:) = [186 384 240];
Small(12,:) = [179 324 245];
Small(13,:) = [118 637 85];
Small(14,:) = [391 506 76];
Small(15,:) = [176 617 68];
Small(16,:) = [137 635 48];
Small(17,:) = [308 367 508];
Small(18,:) = [342 412 376];
Small(19,:) = [331 381 364];
Small(20,:) = [367 535 333];
Small(21,:) = [185 507 410];
Small(22,:) = [82	624	274];

% line that bisects Mi3 and Ca

MeanBorder6_7 = mean([MeanM3C;MeanCa]);
MeanBorder5_6 = mean([MeanM2C; MeanM3C]);
MeanBorder4_5 = mean([MCell;MeanM2C]);

MeanBorder6_7_projcet = Projection(MeanBorder6_7,p1,p3); 
MeanBorder5_6_projcet = Projection(MeanBorder5_6,p1,p3); 
MeanBorder4_5_projcet = Projection(MeanBorder4_5,p1,p3); 

%generate YZ image volume
for i = 1:imWidth
    for j = 1:imNumber
        SideImage(:,j,i) = FinalImage(:,i,j);
    end
end

% calculate border lines
index =2;
x = 1:22:imNumber;
y(1) =  MeanBorder6_7_projcet(2);
for i = 1:length(x)
    y(index) = y(index-1)-1.26; % 370/292 nm
    index = index+1;
end

% plot planes ever 10um; 10*1000/292
figure()
index =1;
for i =  1:35:imWidth
    subplot(2,8,index);
    imshow(medfilt2((SideImage(:,:,i))));
    daspect([6.48,1,1]);
    hold on;
    %scatter(x,y(1:80),'k.');
    index = index+1;
    title(i*0.292);
    scatter(MCell(3), MCell(2),500,'k', 'p', 'MarkerFaceColor','k');
    scatter(MeanM2C(3),MeanM2C(2),300,'r', 'p', 'MarkerFaceColor','r', 'MarkerEdgeColor','k');
    scatter(MeanM3C(3),MeanM3C(2),300,'b','p', 'MarkerFaceColor','b', 'MarkerEdgeColor','k');
    scatter(MeanCa(3), MeanCa(2),300, 'g','p', 'MarkerFaceColor','g', 'MarkerEdgeColor','k');
    line([1500 1500], [100, 168.49],'Color','w', 'LineWidth', 4) % 20 um line
   % scatter(Small(:,3),Small(:,2),100, 'm','o','MarkerFaceColor','m', 'MarkerEdgeColor','k');
end


% plot line across planes with abducens nerve

figure();
subplot(1,2,1);
imshow(medfilt2((SideImage(:,:,213))));daspect([6.48,1,1]);
hold on;
scatter(MCell(3), MCell(2),500,'k', 'p', 'MarkerFaceColor','k');
scatter(MeanM2C(3),MeanM2C(2),300,'r', 'p', 'MarkerFaceColor','r', 'MarkerEdgeColor','k');
scatter(MeanM3C(3),MeanM3C(2),300,'b','p', 'MarkerFaceColor','b', 'MarkerEdgeColor','k');
scatter(MeanCa(3), MeanCa(2),300, 'g','p', 'MarkerFaceColor','g', 'MarkerEdgeColor','k');
scatter(Small(:,3),Small(:,2),100, 'm','o','MarkerFaceColor','m', 'MarkerEdgeColor','k');
scatter(x,y(1:80),'k.');
line([1500 1500], [100, 168.49],'Color','w', 'LineWidth', 4) % 20 um line
line([1111 1111], [903 0 ], 'Color','w', 'LineWidth', 4, 'LineStyle', '--');
title('Rostral abducens nerve');

subplot(1,2,2)
imshow(medfilt2((SideImage(:,:,231))));daspect([6.48,1,1]);
hold on;
scatter(MCell(3), MCell(2),500,'k', 'p', 'MarkerFaceColor','k');
scatter(MeanM2C(3),MeanM2C(2),300,'r', 'p', 'MarkerFaceColor','r', 'MarkerEdgeColor','k');
scatter(MeanM3C(3),MeanM3C(2),300,'b','p', 'MarkerFaceColor','b', 'MarkerEdgeColor','k');
scatter(MeanCa(3), MeanCa(2),300, 'g','p', 'MarkerFaceColor','g', 'MarkerEdgeColor','k');
scatter(Small(:,3),Small(:,2),100, 'm','o','MarkerFaceColor','m', 'MarkerEdgeColor','k');
scatter(x,y(1:80),'k.');
line([1500 1500], [100, 168.49],'Color','w', 'LineWidth', 4) % 20 um line
line([1111 1111], [0 903], 'Color','w', 'LineWidth', 4, 'LineStyle', '--');
title('Caudal abducens nerve');






