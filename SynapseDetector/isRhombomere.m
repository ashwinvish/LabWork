function [inBoundary] = isRhombomere(cellID)
% isRhombomere deteremines which rhombomere the root node of cellID
% lies in.

% path to relavent folder with the skeletons.

if ismac
    addpath(genpath('/Users/ashwin/Documents/LabWork'));
    fname  = '/Users/ashwin/Documents/LowEMtoHighEM/SWC_all/consensus-20180920/swc/';
    load('/Users/ashwin/Documents/RefBrains/MaskDatabase.mat');
    imagesFile = '/Users/ashwin/Documents/RefBrains/SpinalBackfills.tif';
    
else
    addpath(genpath('/usr/people/ashwinv/seungmount/research/Ashwin/Scripts'));
    fname =  '/usr/people/ashwinv/seungmount/research/Ashwin/Z-Brain/LowEMtoHighEM/SWC_all/consensus-20180920/swc/'
    load('/usr/people/ashwinv/seungmount/research/Ashwin/Z-Brain/RefBrains/MaskDatabase.mat');
    imagesFile = '/usr/people/ashwinv/seungmount/research/Ashwin/Z-Brain/RefBrains/SpinalBackfills.tif';
end


% need these params to set the global scale
imagesInfo = imfinfo(imagesFile);
width = imagesInfo(1).Width;
height = imagesInfo(1).Height;
imagesRef = imref3d([width,height,length(imagesInfo)],0.798,0.798,2); % set the correct dimensions for the images using the ref object


IdsHighlight = [221:225];        % r3:r7 ;

for i = 1:length(IdsHighlight)
    maskImageHighlight(i).Ids = IdsHighlight(i);
    maskImageHighlight(i).name = MaskDatabaseNames(IdsHighlight(i));
    maskImageHighlight(i).image = reshape(full(MaskDatabase(:,IdsHighlight(i))),height,width,length(imagesInfo));
end

%  load and transfrom skeleton to z-brain atlas space

%swc_zbrian = zeros(1,size(cellID,1));
rootNodePlane = zeros(size(cellID,1),3);
for i = 1:length(cellID)
    swc_zbrian{i} = SwctoZbrian(cellID(i));
    rootNodePlane(i,:) = [swc_zbrian{i}{1}.X(1),swc_zbrian{i}{1}.Y(1),swc_zbrian{i}{1}.Z(1)]; % get root node coords in Z-brain space.
end

% determine which rhombomere the rootnode lies in.
for j = 1:length(cellID)
    for i = 1:length(IdsHighlight)
        invertedPlaneinMicrons = imagesRef.ImageExtentInWorldZ-rootNodePlane(j,3);
        invertedPlaneinVoxel = round(invertedPlaneinMicrons/imagesRef.PixelExtentInWorldZ);
        [B,L] = bwboundaries(maskImageHighlight(i).image(:,:,invertedPlaneinVoxel));
        inveratedPlaneinMicrons = rootNodePlane(j,3);
        temp = cell2mat(B);
        if i ==1
            boundaries(j).r3(:,1) =  temp(:,1).*[imagesRef.PixelExtentInWorldX];
            boundaries(j).r3(:,2) = temp(:,2).*[imagesRef.PixelExtentInWorldY];
            boundaries(j).r3(:,3) = invertedPlaneinMicrons*ones(size(temp,1),1);
            inBoundary(j).r3 = inpolygon(rootNodePlane(j,1),rootNodePlane(j,2),boundaries(j).r3(:,2), boundaries(j).r3(:,1));
        elseif i ==2
            boundaries(j).r4(:,1) =  temp(:,1).*[imagesRef.PixelExtentInWorldX];
            boundaries(j).r4(:,2) = temp(:,2).*[imagesRef.PixelExtentInWorldY];
            boundaries(j).r4(:,3) = invertedPlaneinMicrons*ones(size(temp,1),1);
            inBoundary(j).r4 = inpolygon(rootNodePlane(j,1),rootNodePlane(j,2),boundaries(j).r4(:,2), boundaries(j).r4(:,1));

        elseif i ==3
            boundaries(j).r5(:,1) =  temp(:,1).*[imagesRef.PixelExtentInWorldX];
            boundaries(j).r5(:,2) = temp(:,2).*[imagesRef.PixelExtentInWorldY];
            boundaries(j).r5(:,3) = invertedPlaneinMicrons*ones(size(temp,1),1);
            inBoundary(j).r5 = inpolygon(rootNodePlane(j,1),rootNodePlane(j,2),boundaries(j).r5(:,2), boundaries(j).r5(:,1));

        elseif i ==4
            boundaries(j).r6(:,1) =  temp(:,1).*[imagesRef.PixelExtentInWorldX];
            boundaries(j).r6(:,2) = temp(:,2).*[imagesRef.PixelExtentInWorldY];
            boundaries(j).r6(:,3) = invertedPlaneinMicrons*ones(size(temp,1),1);
            inBoundary(j).r6 = inpolygon(rootNodePlane(j,1),rootNodePlane(j,2),boundaries(j).r6(:,2), boundaries(j).r6(:,1));

        elseif i ==5
            boundaries(j).r7(:,1) =  temp(:,1).*[imagesRef.PixelExtentInWorldX];
            boundaries(j).r7(:,2) = temp(:,2).*[imagesRef.PixelExtentInWorldY];
            boundaries(j).r7(:,3) = invertedPlaneinMicrons*ones(size(temp,1),1);
            inBoundary(j).r7 = inpolygon(rootNodePlane(j,1),rootNodePlane(j,2),boundaries(j).r7(:,2), boundaries(j).r7(:,1));
        end  
    end

end
end

