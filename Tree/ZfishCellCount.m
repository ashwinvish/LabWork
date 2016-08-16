%% Load images
%ImageFname = '/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/Experiements/10122012-1/ZFishLowResAlined1285.tif'; % file name for the stack used
ImageFname = '/Users/admin/Documents/ZFishLowResAlined1285.tif';
info = imfinfo(ImageFname);
for i = 1:numel(info)
    %ZFishImage(:,:,i) = imread(ImageFname,i,'PixelRegion',{[4, 754], [34,427]});
    ZFishImage(:,:,i) = imread(ImageFname,i,'PixelRegion',{[0, 903], [0, 511]});
end
%% local cellbody location

clear ZfishCellCounts;
clear ZfishSmallX;
clear ZfishSmallY;
clear ZfishSmallZ;

fileID = fopen('Results.txt');
ZfishCellCounts = textscan(fileID, '%d %d %d %d %f','Delimiter','\t');
fclose(fileID);

% display only those that are in the High-res volume
index =1;
% for i = 1:length(ZfishCellCounts{1,3})
%     if ((34> ZfishCellCounts{1,3}(i) <= 427) &&  (4> ZfishCellCounts{1,4}(i) <= 732)) ==1
%         ZfishSmallX(index) = ZfishCellCounts{1,3}(i);
%         ZfishSmallY(index) = ZfishCellCounts{1,4}(i);
%         ZfishSmallZ(index) = -1*ZfishCellCounts{1,2}(i);
%         index = index+1;
%     end
% end

for i = 1:length(ZfishCellCounts{1,3})
    if ZfishCellCounts{1,3}(i)> 34 &&  ZfishCellCounts{1,3}(i) <= 427
        if  ZfishCellCounts{1,4}(i) >4 && ZfishCellCounts{1,4}(i)<= 732
            ZfishSmallX(index) = ZfishCellCounts{1,3}(i);
            ZfishSmallY(index) = ZfishCellCounts{1,4}(i);
            ZfishSmallZ(index) = -1*ZfishCellCounts{1,2}(i);
            index = index+1;
        end
    end
end

%MarkerColors = repmat([0.9,0.9,0.9],length(ZfishSmallX),1);
MarkerColors = repmat([0.1,0.1,0.1],length(ZfishSmallX),1);
hplot = scatter3( ZfishSmallX,ZfishSmallY,ZfishSmallZ ,10,MarkerColors,'o','filled','MarkerEdgeColor','k');
daspect([1,1,9]);
axis vis3d;
set(gca,'XDir','reverse');
%  pause(3);
%  hMarker = hplot.MarkerHandle;
%  hMarker.FaceColorData = uint8(255*[0.9;0.9;0.9;0.3]);
hold on;

% Anaotmical cells of interest
MCell = [210,	828,	458];
M2C1 = 	[117,	736,	533];
M2C2 = 	[150,	748,	398];
M3C1 = 	[189,	653,	47];
M3C2 = 	[196,	670,	115];
CaD =   [80,    583,    1];
CaV =   [86,    604,    23];

scatter3(MCell(1),MCell(2),-1*MCell(3),500,'k', 'p', 'MarkerFaceColor','k');
text( 520,MCell(2),-1, 'r4', 'FontName', 'Arial', 'FontSize', 40, 'Color','k');
%text(-10, MCell(2), -1, 'Mauthner cell','FontName', 'Arial', 'FontSize', 40, 'Color','k');

scatter3(M2C1(1),M2C1(2),-1*M2C1(3),300,'r', 'p', 'MarkerFaceColor','r', 'MarkerEdgeColor','k');
scatter3(M2C2(1),M2C2(2),-1*M2C2(3),300,'r', 'p', 'MarkerFaceColor','r', 'MarkerEdgeColor','k');
MeanM2C = mean([M2C1; M2C2]);
text(520, MeanM2C(2), -1,'r5', 'FontName', 'Arial', 'FontSize', 40, 'Color','k');
%text(-10, MeanM2C(2), -1,'Mi2', 'FontName', 'Arial', 'FontSize', 40, 'Color','r');

scatter3(M3C1(1),M3C2(2),-1*M2C2(3),300,'b', 'p', 'MarkerFaceColor','b', 'MarkerEdgeColor','k');
scatter3(M3C2(1),M3C2(2),-1*M3C2(3),300,'b','p', 'MarkerFaceColor','b', 'MarkerEdgeColor','k');
MeanM3C = mean([M3C1; M3C2]);
text(520, MeanM3C(2), -1, 'r6', 'FontName', 'Arial', 'FontSize', 40, 'Color','k');
%text(-10, MeanM3C(2), -1, 'Mi3', 'FontName', 'Arial', 'FontSize', 40, 'Color','b');

scatter3(CaD(1), CaD(2), -1*CaD(3), 300, 'g','p', 'MarkerFaceColor','g', 'MarkerEdgeColor','k');
scatter3(CaV(1), CaV(2), -1*CaV(3), 300, 'g','p', 'MarkerFaceColor','g', 'MarkerEdgeColor','k');
MeanCa  = mean([CaD; CaV]);
text(520, MeanCa(2), -1, 'r7', 'FontName', 'Arial', 'FontSize', 40, 'Color','k');
%text(-10, MeanCa(2), -1, 'Ca', 'FontName', 'Arial', 'FontSize', 40, 'Color','g');

text(520,412,-1,'M1','FontName', 'Arial', 'FontSize', 40, 'Color','k');
text(520,270,-1,'M2','FontName', 'Arial', 'FontSize', 40, 'Color','k');
text(520,109,-1,'M3','FontName', 'Arial', 'FontSize', 40, 'Color','k');

set(gca, 'color','none');
box on;
set(gca, 'BoxStyle','full', 'XTick',[],'YTick',[], 'ZTick',[]);

% RC distance between reticulospinal cells

% creating a flat surface of same size as image
%colormap gray;
%[X,Y] = meshgrid(linspace(0, 393, 393), linspace(0,750,750)); % image dimensions are 903x511
[X,Y] = meshgrid(linspace(0, 511, 511), linspace(0,903,903)); % image dimensions are 903x511
imageplane = 458; % 1 to 1285
Z = -(imageplane)*ones(size(X));
hsurf = surface(X,Y,Z,histeq(ZFishImage(:,:,imageplane),120),'FaceColor','texturemap','EdgeColor','none')
alpha(hsurf, 0.7);
% alpha color;
% alpha scaled;
set(gca, 'BoxStyle','full')
view(30,60);
axis off;

%% Add VPNI cells
hold on;
SmallAlx = [4,5,6,13,16,22];
SmallTrans = [7,21];
SmallDbx = [2,3,8,9,10,11,12,15];
SmallBarhl = [1,14,17,18,19,20];

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

% for i = 1:length(Small)
%     if ismember(i,SmallAlx) ==1
%         hcell = scatter3(Small(i,1),Small(i,2),-1*Small(i,3),250, calx,'filled','MarkerEdgeColor','k');
%         uistack(hcell,'top');
%         pause(1);
%        MarkerColorMap(:,i) = calx';
%     elseif ismember(i, SmallTrans) ==1
%        hcell =  scatter3(Small(i,1),Small(i,2),-1*Small(i,3),250, ctrans,'filled','MarkerEdgeColor','k');
%         uistack(hcell,'top');
%          pause(1);
%         MarkerColorMap(:,i) = ctrans';
%     elseif ismember( i, SmallDbx) ==1
%         hcell = scatter3(Small(i,1),Small(i,2),-1*Small(i,3),250, cdbx,'filled','MarkerEdgeColor','k');
%          uistack(hcell,'top');
%           pause(1);
%         MarkerColorMap(:,i) = cdbx';
%     else
%         hcell = scatter3(Small(i,1),Small(i,2),-1*Small(i,3),250, cbarhl,'filled','MarkerEdgeColor','k');
%          uistack(hcell,'top');
%           pause(1);
%         MarkerColorMap(:,i) = cbarhl';
%     end
% end
 CT = cbrewer('div','Spectral',size(CellSoma,1));
% %CT = parula(22);
 [y,I] = sort(rho);
 [y,I] =  sort(log2(tau));
 for i = 1:numel(cellIDs)
     scatter3(Small(I(i),1), Small(I(i),2), Small(I(i),3), 250, CT(i,:), 'filled', 'MarkerEdgeColor','k');
 end
%colorbar;
daspect([1,1,9]);
set(gca,'XDir','reverse');
colormap(CT);
%MarkerColorMap = [MarkerColorMap;repmat(1,1,22)];

%% Cell counts to get peaks
figure();
scatter3( ZfishSmallX,ZfishSmallY,ZfishSmallZ ,50,MarkerColors,'filled','MarkerEdgeColor','k');
daspect([1,1,9]);
axis vis3d;
set(gca,'XDir','reverse');
view (0,90)
hold on;
% creating a flat surface of same size as image
colormap gray;
[X,Y] = meshgrid(linspace(0, 393, 393), linspace(0,750,750)); % image dimensions are 903x511
Z = -(imageplane)*ones(size(X));
hsurf = surface(X,Y,Z,ZFishImage(:,:,imageplane),'FaceColor','texturemap','EdgeColor','none')
alpha(hsurf, 0.7);
set(gca, 'BoxStyle','full');
%axis off;

%%
figure();

% h = histogram(ZfishSmallX(find(-1*ZfishSmallZ<imageplane) && find(ZfishSmallY> 200)),'BinWidth',4,'Visible','off');
%  CellCount = [0,h.Values];
%  CellXAxis = h.BinEdges;
%  hold on
% plot(0.293*CellXAxis, smooth(CellCount),'-k', 'LineWidth',2);
% box off;
% xlabel('Distance (\mum)', 'FontName', 'Arial', 'FontSize', 40);
% ylabel('Number of cells', 'FontName', 'Arial', 'FontSize', 40);
% set(gca, 'XDir','reverse','FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
index = 1;
for i = 1:length(ZfishSmallX)
    if ((ZfishSmallY (i) > 150) &&  (-1*ZfishSmallZ (i)<imageplane) )==1
        XStripe(1, index) = ZfishSmallX(i);
        index = index+1;
    end
end
h = histogram(XStripe,'BinWidth',4,'Visible','off');
CellCount = [0,h.Values];
CellXAxis = h.BinEdges;

hold on;

% cell soma colored by type

% for i = 1:length(Small)
%     if ismember(i,SmallAlx) ==1
%         [a,h1,h2] = plotyy(0.293*CellXAxis, smooth(CellCount),0.293* Small(i,1),45*Small(i,3) /1000,'line','scatter');
%         set(h1,'LineStyle','none');
%         set(h2,'MarkerFaceColor',calx,'MarkerEdgeColor','k', 'SizeData',2000);
%         set(a(2), 'XDir','reverse', 'XTick',[], 'YLim',[0 25], 'YTick',[], 'YDir','reverse');
%         hold on;
%     elseif ismember(i, SmallTrans) ==1
%         [a,h1,h2] = plotyy(0.293*CellXAxis, smooth(CellCount),0.293*Small(i,1),45*Small(i,3) /1000,'line','scatter');
%         set(h1,'LineStyle','none');
%         set(h2,'MarkerFaceColor',ctrans,'MarkerEdgeColor','k', 'SizeData',2000);
%         set(a(2), 'XDir','reverse', 'XTick',[], 'YLim',[0 25],'YTick',[], 'YDir','reverse');
%         hold on;
%     elseif ismember( i, SmallDbx) ==1
%         [a,h1,h2] = plotyy(0.293*CellXAxis, smooth(CellCount),0.293*Small(i,1),45*Small(i,3) /1000,'line','scatter');
%         set(h1,'LineStyle','none');
%         set(h2,'MarkerFaceColor',cdbx,'MarkerEdgeColor','k', 'SizeData',2000);
%         set(a(2), 'XDir','reverse', 'XTick',[],'YLim',[0 25], 'YTick',[], 'YDir','reverse');
%         hold on;
%     else
%         [a,h1,h2] = plotyy(0.293*CellXAxis, smooth(CellCount),0.293*Small(i,1),45*Small(i,3) /1000,'line','scatter')
%         set(h1,'LineStyle','none');
%         set(h2,'MarkerFaceColor',cbarhl,'MarkerEdgeColor','k', 'SizeData',2000);
%         set(a(2), 'XDir','reverse', 'XTick',[],'YLim',[0 25],'YTick',[], 'YDir','reverse');
%         hold on;
%     end
% end

% cell soma colored by rhos

[y,I] = sort(rho);
for i = 1:numel(cellIDs)
    [a,h1,h2] = plotyy(0.293*CellXAxis, smooth(CellCount),0.293* Small(I(i),1),45*Small(I(i),3) /1000,'line','scatter');
    set(h1,'LineStyle','none');
    set(h2,'MarkerFaceColor',CT(i,:),'MarkerEdgeColor','k', 'SizeData',5000);
    set(a(2), 'XDir','reverse', 'XTick',[],'YLim',[0 25], 'YTick',[], 'YDir','reverse');
    hold on;
end


[pks,loc] = findpeaks(smooth(CellCount),'MinPeakHeight',17);
[a,h1,h2] = plotyy(0.293*CellXAxis, smooth(CellCount),0,0,'line','scatter');
set(h1,'Color','b', 'LineWidth', 4);
set(h2,'MarkerFaceColor','none','MarkerEdgeColor','none');

% insert stripe locations, 0.293*CellXAxis(loc)
plot(a(1), [0.293*CellXAxis(loc(2)),0.293*CellXAxis(loc(2)) ], [0, 25], 'k--', 'LineWidth',4);
plot(a(1), [0.293*CellXAxis(loc(3)),0.293*CellXAxis(loc(3)) ], [0, 25], 'k--','LineWidth',4);
plot(a(1), [0.293*CellXAxis(loc(4)),0.293*CellXAxis(loc(4)) ], [0, 25], 'k--','LineWidth',4);


set(a(1),'XDir','reverse','XLim',[0, 140],'YColor','b','YLim',[0,25],'YTick',[0,5,10,15,20,25],'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2 );
set(a(2), 'XDir','reverse', 'XTick',[],'YDir','reverse', 'YLim',[0, 25],'YTick',[0,5,10,15,20,25],'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
xlabel('Distance (\mum)', 'FontName', 'Arial', 'FontSize', 40);
ylabel(a(1),'Number of cells', 'FontName', 'Arial', 'FontSize', 40);
ylabel(a(2),'Depth (\mum)', 'FontName', 'Arial', 'FontSize', 40);




