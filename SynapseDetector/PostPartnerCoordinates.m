function PostCoordinates = PostPartnerCoordinates(PostPSDSynapseSegID,df)
% *PostPSDCoordinates* is used to obtain the x,y,z coordinates in NG
% framework for the listed PostPartnerSynapseID

PostX = [];
PostY = [];
PostZ = [];

    for i = 1: length(PostPSDSynapseSegID)
         PostX = [PostX; df.postsyn_x(df.psd_segid==PostPSDSynapseSegID(i))];
         PostY = [PostY; df.postsyn_x(df.psd_segid==PostPSDSynapseSegID(i))];
         PostZ = [PostZ; df.postsyn_x(df.psd_segid==PostPSDSynapseSegID(i))];
    end
    PostCoordinates = [PostX,PostY,PostZ];
end
