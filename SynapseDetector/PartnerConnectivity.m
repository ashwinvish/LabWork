function [leadLagAxon,ConnLeadLag,LeadLagGraph,LeadLagGraphMetrics] = PartnerConnectivity(leadLagAxonIDs,allUniqueAxonIDs,df)



for i = 1:numel(leadLagAxonIDs)
        leadLagAxon(i) = InputsByClass(leadLagAxonIDs(i),df);
end

%leadLagAxonOrder = find(ismember(allUniqueAxonIDs,leadLagAxonIDs));
%ConnLeadLag = zeros(size(leadLagAxonOrder,1),size(allUniqueAxonIDs,1));

for i = 1:numel(leadLagAxon)
     ConnLeadLag(i,:) = histc(leadLagAxon(i).Saccadic, allUniqueAxonIDs);
end
%ConnLeadLag = ConnLeadLag(:,leadLagAxonOrder);
%ConnLeadLag = ConnLeadLag.*~eye(size(ConnLeadLag));

% graph fromalism
[m,n,v] = find(ConnLeadLag);                    % m -postSynapse,row, n - column, v - number of synapses,value
LeadLagGraph = digraph(n',m',v');                 % digraph(source, target, weight, numNodes)

index = 1;
for i = 1:numel(leadLagAxon)
    for j = 1:numel(leadLagAxon)
        if ~isempty(shortestpath(LeadLagGraph,i,j))
            [p(index).path,p(index).dist] = shortestpath(LeadLagGraph,i,j);         
            p(index).ij = [i,j];
            p(index).size = size(p(index).path,2);
            index = index+1;
        end
    end
end
LeadLagGraphMetrics = p;
end