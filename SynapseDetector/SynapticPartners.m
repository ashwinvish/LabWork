function [PartnerSeg, PartnerCoordinate] = SynapticPartners(cellID, partnerID, df)
% SynapticPartners finds all the synaptically connected partners 
% cellID is the ID of the cells for which partners are queried
% partnerID is presynapse = 1 , postsynapse = 2
% df is the dataframe
% partners is mx2 returns the synapticpartnerID and the PSD_ID.

if partnerID == 1   
    PartnerSeg              = df.presyn_seg(df.postsyn_seg==cellID);
    PartnerCoordinate       = df.psd_segid(df.postsyn_seg==cellID);
else
    PartnerSeg              = df.postsyn_seg(df.presyn_seg==cellID);
    PartnerCoordinate(:)    = df.psd_segid(df.presyn_seg==cellID);
end
        
end