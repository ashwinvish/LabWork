function [ SynDistance, minimumDistance, Irow, Icol  ] = MinSynapticDistance( PreSynapses, PostSynapses )
    % MINSYNAPTICDISTANCE minimum pairwise synaptic distance between
    % synapses of two cells
    %   PRESYNAPSES is a Nx3 vector with cartesian coordinates of all
    %   presynaptic sites
    %   POSTSYNAPSES is a Mx3 vector with cartesian coordinates of all
    %   postsynaptic sites
    %   SYNDISTANCE is MxN matrix of all pairwise distances
    %   MINIMUMDISTANCE is the minimum of DISTANCE
    %   [IROW ICOL] is the locaion of the MINIMUMDISTACE in DISTANCE
    
%     for k = 1: length(PreSynapses)
%         SynDistance(1:length(PostSynapses),k) = sqrt((PreSynapses(k,1) - PostSynapses(:,1)).^2 + (PreSynapses(k,2) - PostSynapses(:,2)).^2 ...
%             +(PreSynapses(k,3) - PostSynapses(:,3)).^2);
%     end
    
    SynDistance = pdist2(PreSynapses,PostSynapses)';
    
    [minimumDistance,I] = min(SynDistance(:));
    [Irow Icol] =  ind2sub(size(SynDistance),I);
    
end

