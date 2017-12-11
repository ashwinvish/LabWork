function[transPostSynapticPoints,transSomata] = translateTree(postSynapticPoints,currentOrigin,referenceOrigin)
delta = referenceOrigin - currentOrigin;
%delta = referenceOrigin(1,1,1)- postSynapticPoints(:,:,:);
transPostSynapticPoints = zeros(4, length(postSynapticPoints));
for i =  1:length(postSynapticPoints)
    transPostSynapticPoints(:,i) = [1,0,0,delta(1,1); 0,1,0,delta(1,2); 0,0,1,delta(1,3);0,0,0,0]*[postSynapticPoints(i,1);...
        postSynapticPoints(i,2);postSynapticPoints(i,3);1];
    transSomata = [1,0,0,delta(1,1); 0,1,0,delta(1,2); 0,0,1,delta(1,3);0,0,0,0]*[currentOrigin(1,1);...
        currentOrigin(1,2);currentOrigin(1,3);1];
end
transPostSynapticPoints = transPostSynapticPoints';
transPostSynapticPoints = transPostSynapticPoints(:,1:3);

transSomata = transSomata';
transSomata = transSomata(:,1:3);

end