function point_cloud = array2ind(array,segm_id)

    if ~exist('segm_id','var') || isempty(segm_id)
        [point_cloud(:,1), point_cloud(:,2), point_cloud(:,3)] = ind2sub(size(array),find(array));
    else
        [point_cloud(:,1), point_cloud(:,2), point_cloud(:,3)] = ind2sub(size(array),find(array==segm_id));
    end
    
end