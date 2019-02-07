function  transform_SWC_AV(CELL_Id, neuronColor, somataColor, IdsHighlight,  displayRefBrain)

if ismac
    addpath(genpath('/Users/ashwin/Documents/LabWork'));
    fname  = '/Users/ashwin/Documents/LowEMtoHighEM/SWC_all/consensus-20180920/swc/';
    
else
    addpath(genpath('/usr/people/ashwinv/seungmount/research/Ashwin/Scripts'));
    fname =  '/usr/people/ashwinv/seungmount/research/Ashwin/Z-Brain/LowEMtoHighEM/SWC_all/consensus-20180920/swc/'
end

%% Read reference brain or any brain from the ref atlas
% NOTE: Original files needs to be rotated ccw 90 and the zaxis needs to be
% reversed in FIJI only.


if ismac
    imagesFile = '/Users/ashwin/Documents/RefBrains/SpinalBackfills.tif';
else
    imagesFile = '/usr/people/ashwinv/seungmount/research/Ashwin/Z-Brain/RefBrains/SpinalBackfills.tif';
end

imagesInfo = imfinfo(imagesFile);
width = imagesInfo(1).Width;
height = imagesInfo(1).Height;
%images = zeros(height, width,length(imagesInfo),'uint8');
images = zeros(height, width,length(imagesInfo));
TifLink = Tiff(imagesFile,'r');

for i = 1:length(imagesInfo)
    TifLink.setDirectory(i);
    images(:,:,i) = TifLink.read();
end
warning('off','last');
TifLink.close();
imagesRef = imref3d(size(images),0.798,0.798,2); % set the correct dimensions for the images using the ref object

%% Mask Image

disp('Load Mask database');         % mask database can be obtained from the Z-Brain atlas.

if ismac
    load('/Users/ashwin/Documents/RefBrains/MaskDatabase.mat');
else
    load('/usr/people/ashwinv/seungmount/research/Ashwin/Z-Brain/RefBrains/MaskDatabase.mat');
end

IdsBase = [114,221:225]; %          r3-7; Keep this always ON
for i = 1:length(IdsBase)
    maskImage(i).Ids = IdsBase(i);
    maskImage(i).name = MaskDatabaseNames(IdsBase(i));
    maskImage(i).image = reshape(full(MaskDatabase(:,IdsBase(i))),height,width,length(imagesInfo));
end

% other masks that need to be displayed!

if IdsHighlight
    
    %IdsHighlight = [238,235,186];        % Vestibular Clusters
    %Ids =  [Ids,184:190];                % Reticulospinal Clusters
    %Ids =  [Ids,135,186, 246, 250];      % Gad1b-s2, Gly2-s2, VglutS3, vglut2-strip4
    %Ids =  [93,131,130,134];             % Cerebellum
    
    for i = 1:length(IdsHighlight)
        maskImageHighlight(i).Ids = IdsHighlight(i);
        maskImageHighlight(i).name = MaskDatabaseNames(IdsHighlight(i));
        maskImageHighlight(i).image = reshape(full(MaskDatabase(:,IdsHighlight(i))),height,width,length(imagesInfo));
    end
end


%% get swc coordinate
swc = cell(1,size(CELL_Id,2));
swc_new = cell(1,size(CELL_Id,2));
for i = 1:length(CELL_Id)
    filename = sprintf('%d_reRoot_reSample_5000.swc',CELL_Id(i));
    if exist(fullfile(fname,filename))
        fID = fopen(fileName);
        swc{i} = textscan(fID, '%f %f %f %f %f %f %f','HeaderLines',6,'CollectOutput',true);
        fclose(fID);
    elseif exist(fullfile(fname,sprintf('%d.swc',CELL_Id(i))))
        fID = fopen(fileName);
        swc{i} = textscan(fID, '%f %f %f %f %f %f %f','HeaderLines',0,'CollectOutput',true);
        fclose(fID);
    else
        continue;
    end
    clear coord;
    coord = swc{i}(:, 3:5) ;
    
    % compute offset as data in NG is offset from origin
    offset = [920,752,16400] .* [80, 80, 45]; % offset at mip4 number can be obtained from NG .info file
    
    
    coord(:,1) = coord(:,1) - offset(1);
    coord(:,2) = coord(:,2) - offset(2);
    coord(:,3) = coord(:,3) - offset(3);
    
    % adjust voxel size to LM space
    %voxel --> micron
    
    coord = coord ./ 1000;
    coord(:,1) = coord(:,1)/0.798;
    coord(:,2) = coord(:,2)/0.798;
    coord(:,3) = coord(:,3)/2;
    
    % rotation NG data is 90 to the righ as compared to the z-brain atlas
    % Create rotation matrix
    theta = pi/2; % to rotate 90 counterclockwise
    R = [cos(theta) -sin(theta) 0 0; sin(theta) cos(theta) 0 0; ...
        0 0 1 0; 0 0 0 1];
    rotateTForm = affine3d(R);
    coord = transformPointsForward(rotateTForm, coord);
    
    % translate the image after it is rotate to compensate for the rotation
    % about the top left corner
    coord(:,2) = coord(:,2) + 436;
    
    % perform transformation
    if ismac
        load('/Users/ashwin/Documents/LowEMtoHighEM/tformRough-LowEMtoHighEM-set2.mat');
    else
        load('/usr/people/ashwinv/seungmount/research/Ashwin/Z-Brain/LowEMtoHighEM/tformRough-LowEMtoHighEM-set2.mat');
    end
    
    % transform from voxels to microns
    coord(:,1) = coord(:,1) * 0.798 ;
    coord(:,2) = coord(:,2) * 0.798 ;
    coord(:,3) = coord(:,3) * 2 ;
    
    % transform from NG space to Atlas space using a precomputed
    % transfromation matrix.
    
    coord_transformed = transformPointsForward(tformRough, coord);
    swc_new{i}(:,3:5) = coord_transformed;
    % dlmwrite('76748_LM.swc',  swc_new, ' ');
    
    % get mask plane and image (in mircron space, atlas space) through which the root node traverses.
    
    n=1;
    rootNodePlane(i) = round(min(coord_transformed(n,3))); % in micron space
    rootNodeImage(:,:,i) = images(:,:,round(rootNodePlane(i)/imagesRef.PixelExtentInWorldZ));
end
%% construct a grided volume

% deteremine the mean plane through which the root nodes traverses
meanRootNodePlane = mean(rootNodePlane);

% sample above and below root node plane
meanRootNodePlane = meanRootNodePlane;

[X,Y] = meshgrid(linspace(0,ceil(imagesRef.ImageExtentInWorldX), ceil(imagesRef.ImageExtentInWorldX)),...
    linspace(0,ceil(imagesRef.ImageExtentInWorldY),ceil(imagesRef.ImageExtentInWorldY)));

% determine the mean plane image

%meanImage = mean(rootNodeImage,3);
meanImage = images(:,:,round(meanRootNodePlane/imagesRef.PixelExtentInWorldZ));

Z = round(meanRootNodePlane)*ones(size(X));

% resize the mean image to micron space
img = imresize(meanImage,[imagesRef.ImageExtentInWorldY,imagesRef.ImageExtentInWorldX]);

% construct the mean plane as a surface

if displayRefBrain
    hsurf1 = surface(X,Y,Z,imcomplement(medfilt2(imadjust(img/max(img(:))))),'FaceColor','flat','EdgeColor','none');
    %hold on;
    %hsurf2 = surface(X,Y,Z,img/max(img(:)),'FaceColor','flat','EdgeColor','none');
    colormap gray;
    alpha(hsurf1, 0.5);
    %alpha(hsurf2, 0.5);
end

daspect([1,1,1]);
hold on;

% draw boundary of the anatomical region of interst

cols = cbrewer('qual','Accent',length(IdsBase));

%All masks are in voxel dimensions, need to covert them to micron space.
%masks are also inverted, meaning from V-->D.

% Base level maske IDs (used for oveall orientation, e.g. rhombomere
% boundaries)

for i = 1:length(IdsBase)
    invertedPlaneinMicrons = imagesRef.ImageExtentInWorldZ-meanRootNodePlane;
    invertedPlaneinVoxel = round(invertedPlaneinMicrons/imagesRef.PixelExtentInWorldZ);
    [B,L] = bwboundaries(maskImage(i).image(:,:,invertedPlaneinVoxel));
    invertedPlaneinMicrons = meanRootNodePlane;
    %boundarySize = size(B,1);
    if size(B,1)>0
        boundaries = B{1};
        boundaries = boundaries.*[imagesRef.PixelExtentInWorldX,imagesRef.PixelExtentInWorldY];
        boundaries(:,3) = invertedPlaneinMicrons*ones(size(boundaries,1),1);
        %fill3(boundaries(:,2),boundaries(:,1),boundaries(:,3),cols(i,:),'FaceAlpha',0.3,'LineStyle','none','LineWidth',0.5);
        plot3(boundaries(:,2),boundaries(:,1),boundaries(:,3),'Color',[0.8,0.8,0.8],'LineWidth',2);
        clear boundaries;
    end
    if size(B,1)>1
        boundaries = B{2};
        boundaries = boundaries.*[imagesRef.PixelExtentInWorldX,imagesRef.PixelExtentInWorldY];
        boundaries(:,3) = invertedPlaneinMicrons*ones(size(boundaries,1),1);
        % fill3(boundaries(:,2),boundaries(:,1),boundaries(:,3),cols(i,:),'FaceAlpha',0.3,'LineStyle','none','LineWidth',0.5);
        plot3(boundaries(:,2),boundaries(:,1),boundaries(:,3),'Color',[0.8,0.8,0.8],'LineWidth',2);
    end
end

% IDs to be highlighted (e.g. Tangential Vestinular Nucleus)

if IdsHighlight
    index = 0.1;
    for i = 1:length(IdsHighlight)
        invertedPlaneinMicrons = imagesRef.ImageExtentInWorldZ-meanRootNodePlane;
        invertedPlaneinVoxel = round(invertedPlaneinMicrons/imagesRef.PixelExtentInWorldZ);
        [B,L] = bwboundaries(maskImageHighlight(i).image(:,:,invertedPlaneinVoxel));
        invertedPlaneinMicrons = meanRootNodePlane;
        %boundarySize = size(B,1);
        if size(B,1)>0
            boundaries = B{1};
            boundaries = boundaries.*[imagesRef.PixelExtentInWorldX,imagesRef.PixelExtentInWorldY];
            boundaries(:,3) = invertedPlaneinMicrons*ones(size(boundaries,1),1);
            %fill3(boundaries(:,2),boundaries(:,1),boundaries(:,3),cols(i,:),'FaceAlpha',0.3,'LineStyle','none','LineWidth',0.5);
            plot3(boundaries(:,2),boundaries(:,1),boundaries(:,3),'Color',cols(i,:),'LineWidth',2);
            clear boundaries;
        end
        if size(B,1)>1
            boundaries = B{2};
            boundaries = boundaries.*[imagesRef.PixelExtentInWorldX,imagesRef.PixelExtentInWorldY];
            boundaries(:,3) = invertedPlaneinMicrons*ones(size(boundaries,1),1);
            % fill3(boundaries(:,2),boundaries(:,1),boundaries(:,3),cols(i,:),'FaceAlpha',0.3,'LineStyle','none','LineWidth',0.5);
            plot3(boundaries(:,2),boundaries(:,1),boundaries(:,3),'Color',cols(i,:),'LineWidth',2);
            annotation('textbox',[0.42,0.4-index,0.5,0],'EdgeColor','none','String',maskImageHighlight(i).name,'Color',cols(i,:));
        end
        index = index+0.01;
    end
end

set(gca, 'BoxStyle','full','YDir','reverse','ZDir','reverse');
set(gca, 'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[]);
set(gca,'ZLim',[0,imagesRef.ImageExtentInWorldZ], 'XLim',[0,imagesRef.ImageExtentInWorldX], 'YLim', [0,imagesRef.ImageExtentInWorldY]);

%% plot cells in CELL_Id

% for i = 1:length(CELL_Id)
%     filename = sprintf('%d_reRoot_reSample_5000.swc',CELL_Id(i));
%     if exist(fullfile(fname,filename))
%         plot3(swc_new{i}(:,3), swc_new{i}(:,4), swc_new{i}(:,5),'MarkerFaceColor',neuronColor,'MarkerEdgeColor','none'); % plot neuron
%         %[a_new,b_new] = min(swc_new{i}(:,5));
%         %[a,b] = min(swc{i}(:,5));
%         b=1;
%         scatter3(swc_new{i}(b,3), swc_new{i}(b,4), swc_new{i}(b,5),50,'MarkerFaceColor',somataColor,'MarkerEdgeColor','w'); % plot somata
%     elseif exist(fullfile(fname,sprintf('%d.swc',CELL_Id(i))))
%         scatter3(swc_new{i}(:,3), swc_new{i}(:,4), swc_new{i}(:,5),1,'MarkerFaceColor',neuronColor,'MarkerEdgeColor','none'); % plot neuron
%         %[a_new,b_new] = min(swc_new{i}(:,5));
%         %[a,b] = min(swc{i}(:,5));
%         b=1;
%         scatter3(swc_new{i}(b,3), swc_new{i}(b,4), swc_new{i}(b,5),50,'MarkerFaceColor',somataColor,'MarkerEdgeColor','w'); % plot somata
%         continue;
%     end
% end

for i = 1:length(CELL_Id)
    filename = sprintf('%d_reRoot_reSample_5000.swc',CELL_Id(i));
    if exist(fullfile(fname,filename))
        tree = load_tree(fullfile(fname,filename));
        [I,J] = ind2sub(size(tree.dA),find(tree.dA));
        line([swc_new{i}(J,3) swc_new{i}(I,3)]',[swc_new{i}(J,4) swc_new{i}(I,4)]',[swc_new{i}(J,5) swc_new{i}(I,5)]',...
            'Color',neuronColor,'LineWidth',0.15);
        hold on;
        scatter3(swc_new{i}(1,3), swc_new{i}(1,4), swc_new{i}(1,5),35,'MarkerFaceColor',somataColor,...
            'MarkerEdgeColor','k');
        clear I;
        clear J;
        clear tree;
    elseif exist(fullfile(fname,sprintf('%d.swc',CELL_Id(i))))
        tree = load_tree(fullfile(fname,sprintf('%d.swc',CELL_Id(i))));
        [I,J] = ind2sub(size(tree.dA),find(tree.dA));
        line([swc_new{i}(J,3) swc_new{i}(I,3)]',[swc_new{i}(J,4) swc_new{i}(I,4)]',[swc_new{i}(J,5) swc_new{i}(I,5)]','Color',neuronColor);
        hold on;
        scatter3(swc_new{i}(1,3), swc_new{i}(1,4), swc_new{i}(1,5),50,'MarkerFaceColor',somataColor,'MarkerEdgeColor','none'); 
        clear I;
        clear J;
        clear tree;
    else
        continue;
    end
end

%draw scale bar 1px = 0.798 us
% line([500,563],[1280,1280],'Color','m','LineWidth',4);
% text(500,1300,'50um','FontName','Arila','FontSize',10);

axis off;
%title(sprintf('%s,%d',name,138-round(mean


end

