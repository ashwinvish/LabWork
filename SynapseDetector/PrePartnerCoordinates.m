function PreCoordinates = PrePartnerCoordinates(PrePartnerSynapseSegID,df)
% *PrePartnerCoordinates* is used to obtain the x,y,z coordinates in NG
% framework for the listed PrePartnerSynapseID

PreX = [];
PreY = [];
PreZ = [];

    for i = 1: length(PrePartnerSynapseSegID)
         PreX = [PreX; df.presyn_x(df.presyn_seg==PrePartnerSynapseSegID(i))];
         PreY = [PreY; df.presyn_y(df.presyn_seg==PrePartnerSynapseSegID(i))];
         PreZ = [PreZ; df.presyn_z(df.presyn_seg==PrePartnerSynapseSegID(i))];
    end
    PreCoordinates = [PreX,PreY,PreZ];
end