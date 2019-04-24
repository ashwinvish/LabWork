function[pc50,vol1,vol2] = cellID2pointCloud(cellID);

if ismac
    addpath(genpath('/Users/ashwin/Documents/LabWork'));
    fname  = '/Users/ashwin/Documents/LowEMtoHighEM/SWC_all/consensus-20180920/swc/';
    
else
    addpath(genpath('/usr/people/ashwinv/seungmount/research/Ashwin/Scripts'));
    fname =  '/usr/people/ashwinv/seungmount/research/Ashwin/Z-Brain/LowEMtoHighEM/SWC_all/consensus-20180920/swc/'
end
meshFileName = sprintf('/Users/ashwin/Desktop/AbducensMeshes/%d.obj',cellID);
fileID = fopen(meshFileName);
obj = textscan(fileID, '%s %f %f %f');
fclose(fileID);
maxVertex = max(find(cell2mat(obj{1,1}) == 'v'));
pc = pointCloud([obj{1,2}(1:maxVertex),obj{1,3}(1:maxVertex),obj{1,4}(1:maxVertex)]); % convert to pointcloud object
gridStep = 0.5;
pc50 = pcdownsample(pc,'gridAverage',gridStep); % downsample pointcloud by 50%

fileName = sprintf('%d_reRoot_reSample_5000.swc',cellID);
fID = fopen(fullfile(fname,fileName));
swc= textscan(fID, '%f %f %f %f %f %f %f','HeaderLines',6,'CollectOutput',true);

PointCloudOrigin = [swc{1,1}(1,3), swc{1,1}(1,4), swc{1,1}(1,5)]; % get origin in Zbrain atlas space.
 radius = 5000;
 [indices,dists] = findNearestNeighbors(pc,PointCloudOrigin,radius);
 somaPC = pointCloud([pc.Location(indices,1),pc.Location(indices,2),pc.Location(indices,3)]);

radius = 4800;
[indices,dists] = findNeighborsInRadius(pc,PointCloudOrigin,radius);
%somaPC = pointCloud([pc.Location(indices,1),pc.Location(indices,2),pc.Location(indices,3)]);

pcshow(somaPC);

%method1 to compute convex hull volume
[model,inLiers,outLiers] = pcfitsphere(somaPC,4800);
[Hull,vol1] = convhull(somaPC.Location(:,1),somaPC.Location(:,2),somaPC.Location(:,3)); % convert to voxel space
%trimesh(Hull,somaPC.Location(:,1),somaPC.Location(:,2),somaPC.Location(:,3),'FaceAlpha',0.1);

% method2 to compute average radius
radii = pdist2(repmat(PointCloudOrigin,size(somaPC.Location,1),1),somaPC.Location);
meanRadius = mean(diag(radii));
vol2 = 4/3 * pi * (meanRadius)^3;
end