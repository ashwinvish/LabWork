function [pointCloudPost] =  TransformPoints(pointCloudPre)

if ismac
    addpath(genpath('/Users/ashwin/Documents/LabWork'));
    fname  = '/Users/ashwin/Documents/LowEMtoHighEM/SWC_all/consensus-20180920/swc/';
    
else
    addpath(genpath('/usr/people/ashwinv/seungmount/research/Ashwin/Scripts'));
    fname =  '/usr/people/ashwinv/seungmount/research/Ashwin/Z-Brain/LowEMtoHighEM/SWC_all/consensus-20180920/swc/'
end
    coord = pointCloudPre;
    
    % multiply with mip0 resolution to get to nm space
    coord = coord .* [5,5,45]; 
    
     % convert the mip0 offset to nm space
    offset = [14720, 12032, 16400] .* [5,5,45]; % offset can be obtained from NG .info file.
    
    % subtract offset to line up with origin
    coord(:,1) = coord(:,1) - offset(1);
    coord(:,2) = coord(:,2) - offset(2);
    coord(:,3) = coord(:,3) - offset(3);
    
    % adjust voxel size to LM space
    %voxel --> micron
    
    coord = coord ./ 1000;
    coord(:,1) = coord(:,1)/0.798;
    coord(:,2) = coord(:,2)/0.798;
    coord(:,3) = coord(:,3)/2;
    
    % rotate the above points to get to z-brian orientation. NG space is
    % rotated 90 to the right as compared to z-brian atlas.
    
    % Create rotation matrix
    theta = pi/2; % to rotate 90 counterclockwise
    R = [cos(theta) -sin(theta) 0 0; sin(theta) cos(theta) 0 0; ...
        0 0 1 0; 0 0 0 1];
    rotateTForm = affine3d(R);
    coord = transformPointsForward(rotateTForm, coord);
    
    % following rotation a translation needs to be applied to bring back to
    % origin. 
    
    coord(:,2) = coord(:,2) + 436;
    
    % perform transformation (transforms were computed in micron space)
    load('tformRough-LowEMtoHighEM-set2.mat')
    
    % transform from voxels to microns.
    coord(:,1) = coord(:,1) * 0.798 ;
    coord(:,2) = coord(:,2) * 0.798 ;
    coord(:,3) = coord(:,3) * 2 ;
    
    % transform from NG space to z-brainAtlas space using a precomputed
    % transfromation matrix that was loaded above.
    
    pointCloudPost = transformPointsForward(tformRough, coord);
   
end

