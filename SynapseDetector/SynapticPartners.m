function [PartnerSeg, PartnerPSD_ID] = SynapticPartners(cellID, partnerType, df)
% SynapticPartners finds all the synaptically connected partners 
% cellID is the ID of the cells for which partners are queried
% partnerType is presynapse = 1 , postsynapse = 2
% df is the dataframe
% partners is mx2 returns the synapticpartnerID and the PSD_ID.

if partnerType == 1   
    PartnerSeg          = df.presyn_segid(df.postsyn_segid==cellID);
    PartnerPSD_ID       = df.psd_segid(df.postsyn_segid==cellID);
else
    PartnerSeg          = df.postsyn_segid(df.presyn_segid==cellID);
    PartnerPSD_ID(:)    = df.psd_segid(df.presyn_segid==cellID);
end
        
end