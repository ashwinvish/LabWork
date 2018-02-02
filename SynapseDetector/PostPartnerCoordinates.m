function PostCoordinates = PostPartnerCoordinates(PostPartnerSynapseSegID,df)
% *PostPartnerCoordinates* is used to obtain the x,y,z coordinates in NG
% framework for the listed PostPartnerSynapseID

PostX = [];
PostY = [];
PostZ = [];

    for i = 1: length(PostPartnerSynapseSegID)
         PostX = [PostX; df.postsyn_x(df.postsyn_seg==PostPartnerSynapseSegID(i))];
         PostY = [PostY; df.postsyn_x(df.postsyn_seg==PostPartnerSynapseSegID(i))];
         PostZ = [PostZ; df.postsyn_x(df.postsyn_seg==PostPartnerSynapseSegID(i))];
    end
    PostCoordinates = [PostX,PostY,PostZ];
end