function[Vol] = cellID2pointCloud(cellID);

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
pc50 = pcdownsample(pc,'random',0.5); % downsample pointcloud by 50%

fileName = sprintf('%d_reRoot_reSample_5000.swc',cellID);
fID = fopen(fullfile(fname,fileName));
swc= textscan(fID, '%f %f %f %f %f %f %f','HeaderLines',6,'CollectOutput',true);

PointCloudOrigin = [swc{1,1}(1,3), swc{1,1}(1,4), swc{1,1}(1,5)]; % get origin in Zbrain atlas space.
radius = 1900;
[indices,dists] = findNearestNeighbors(pc50,PointCloudOrigin,radius);
somaPC = pointCloud([pc50.Location(indices,1),pc50.Location(indices,2),pc50.Location(indices,3)]);
%[model,inLiers,outLiers] = pcfitsphere(somaPC,2600);
[Hull,Vol] = convhull(somaPC.Location(:,1),somaPC.Location(:,2),somaPC.Location(:,3)); % convert to voxel space
%trimesh(Hull,somaPC.Location(:,1),somaPC.Location(:,2),somaPC.Location(:,3))
end