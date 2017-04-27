% function [Test,seg_raw,seg_user1,seg_user2] = compare_validation(file_raw,file_user1,file_user2)
function output = compare_validation(file_raw,file_user1,file_user2)

fprintf('reading the raw file ...\n');
seg_raw=hdf5read(file_raw,'/main');
fprintf('reading the first user file ...\n');
seg_user1=hdf5read(file_user1,'/main');
fprintf('reading the second user file ...\n');
seg_user2=hdf5read(file_user2,'/main');

lst_segid_in_user1=setdiff(unique(seg_user1(:)),0);
lst_segid_in_user2=setdiff(unique(seg_user2(:)),0);

n_segments_in_user1=numel(lst_segid_in_user1);
n_segments_in_user2=numel(lst_segid_in_user2);

lst_supervoxels_in_segments_in_user1=cell(n_segments_in_user1);
lst_supervoxels_in_segments_in_user2=cell(n_segments_in_user2);

lst_supervoxels_in_raw=setdiff(unique(seg_raw(:)),0);
n_supervoxels=numel(lst_supervoxels_in_raw);
size_supervoxels=zeros(1,n_supervoxels);
Test = zeros(1,3);

fprintf('estimating supervoxel sizes and user validation info ...\n');
for i=1:n_supervoxels
    supervoxel_id=lst_supervoxels_in_raw(i);
    size_supervoxels(i)=sum(seg_raw(:)==supervoxel_id);
    sub=find(seg_raw(:)==supervoxel_id,1);

    segid1=seg_user1(sub);
    idx=find(lst_segid_in_user1==segid1);
    lst_supervoxels_in_segments_in_user1{idx}=[lst_supervoxels_in_segments_in_user1{idx} supervoxel_id];

    segid2=seg_user2(sub);
    idx=find(lst_segid_in_user2==segid2);
    lst_supervoxels_in_segments_in_user2{idx}=[lst_supervoxels_in_segments_in_user2{idx} supervoxel_id];
 
    % fprintf('%d %d %d %d\n',supervoxel_id,size_supervoxels(i),segid1,segid2);
end

output=zeros(n_segments_in_user1,n_segments_in_user2);
for i=1:n_segments_in_user1
    lst_overlap_size=zeros(1,n_segments_in_user2);
    for j=1:n_segments_in_user2
        overlap_supervoxels=intersect(lst_supervoxels_in_segments_in_user1{i},lst_supervoxels_in_segments_in_user2{j});
        lst_overlap_size(j)=sum(size_supervoxels(ismember(lst_supervoxels_in_raw,overlap_supervoxels)));
        output(i,j)=sum(size_supervoxels(ismember(lst_supervoxels_in_raw,overlap_supervoxels)));
    end
end

% for i=1:n_segments_in_user1
%     lst_overlap_size=zeros(1,n_segments_in_user2);
%     for j=1:n_segments_in_user2
%         overlap_supervoxels=intersect(lst_supervoxels_in_segments_in_user1{i},lst_supervoxels_in_segments_in_user2{j});
%         lst_overlap_size(j)=sum(size_supervoxels(ismember(lst_supervoxels_in_raw,overlap_supervoxels)));
%     end
%     [lst_overlap_size,IX]=sort(lst_overlap_size,'descend');
% 
%     %fprintf('[% 4d] in user1 correspond user2''s ...\n',lst_segid_in_user1(i));
%     for j=1:n_segments_in_user2
%         if (lst_overlap_size(j)==0)
%             break;
%         end
%         
%         %fprintf(' -[% 4d] by %d\n',lst_segid_in_user2(IX(j)),lst_overlap_size(j));
%         temp = [lst_segid_in_user1(i),lst_segid_in_user2(IX(j)),lst_overlap_size(j)];
%         Test = vertcat(Test,temp);
%         %dlmwrite('test.txt',temp,'delimiter', '\t','newline','unix','-append');
%     end
%     %sprintf(
%     %fprintf('\n');
% end


