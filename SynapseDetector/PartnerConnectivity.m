function [leadLagAxon,ConnLeadLag,LeadLagGraph,LeadLagGraphMetrics] = PartnerConnectivity(leadLagAxonIDs,allUniqueAxonIDs,neuronType,df)


for i = 1:numel(leadLagAxonIDs)
    leadLagAxon(i) = InputsByClass(leadLagAxonIDs(i),df);
end

%leadLagAxonOrder = find(ismember(allUniqueAxonIDs,leadLagAxonIDs));
%ConnLeadLag = zeros(size(leadLagAxonOrder,1),size(allUniqueAxonIDs,1));

if contains('Saccadic',neuronType) == 1
    for i = 1:numel(leadLagAxon)
        ConnLeadLag(i,:) = histc(leadLagAxon(i).Saccadic, allUniqueAxonIDs);
    end
elseif contains('Integrator', neuronType) == 1
    for i = 1:numel(leadLagAxon)
        ConnLeadLag(i,:) = histc(leadLagAxon(i).Integrator, allUniqueAxonIDs);
    end
elseif contains('Contra',neuronType) == 1
    for i = 1:numel(leadLagAxon)
        ConnLeadLag(i,:) = histc(leadLagAxon(i).Contra, allUniqueAxonIDs);
    end
elseif contains('Vestibular',neuronType) == 1
    for i = 1:numel(leadLagAxon)
        ConnLeadLag(i,:) = histc(leadLagAxon(i).Vestibular, allUniqueAxonIDs);
    end
else
    for i = 1:numel(leadLagAxon)
        ConnLeadLag(i,:) = histc(leadLagAxon(i).EverythingElse, allUniqueAxonIDs);
    end
end




%ConnLeadLag = ConnLeadLag(:,leadLagAxonOrder);
%ConnLeadLag = ConnLeadLag.*~eye(size(ConnLeadLag));

% graph fromalism
[m,n,v] = find(ConnLeadLag);                    % m -postSynapse,row, n - column, v - number of synapses,value
LeadLagGraph = digraph(n',m',v');                 % digraph(source, target, weight, numNodes)

index = 1;

for i = 1:numel(leadLagAxon)
    for j = 1:numel(leadLagAxon)
        index = sub2ind([numel(leadLagAxon), numel(leadLagAxon)],i,j);
        if ~isempty(shortestpath(LeadLagGraph,i,j))
            [p(index).path,p(index).dist] = shortestpath(LeadLagGraph,i,j);
            p(index).ij = [i,j];
            p(index).size = size(p(index).path,2);
        else
            p(index).path = [];
            p(index).dist = [];
            p(index).ij = [i,j];
            p(index).size = [];
        end
    end
end
LeadLagGraphMetrics = p;
end