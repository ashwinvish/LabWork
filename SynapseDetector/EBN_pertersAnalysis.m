% EBN Peters rule
clear;

load ABDPutativeSaccadic.mat
load ABDiPutativeSaccadic.mat
load ABDr.mat
load ABDc.mat
load ABDIr.mat
load ABDIc.mat
LoadDataFrame

AllABDPreSynapticLocaions = [vertcat(ABDr.PreSynCoords); vertcat(ABDc.PreSynCoords)];
AllABDiPSDlocs = [vertcat(ABDIr.PSDID); vertcat(ABDIr.PSDID)];
AllABDiPostSynapticLocations = PostPartnerCoordinates(AllABDiPSDlocs(1:100:end),df);


prePartnercoords = [];
postPartnercoords = [];
for i = 1:size(ABDPutativeSaccadic.cellIDs,2)
    [partnerID,partnerPSD] = SynapticPartners(ABDPutativeSaccadic.cellIDs(i),2,df);
    prePartnercoords = [prePartnercoords;PrePartnerCoordinates(partnerPSD,df)];
    postPartnercoords = [postPartnercoords;PostPartnerCoordinates(partnerPSD,df)];
    clear partnerID;
    clear partnerPSD;
end

distMat = pdist2(prePartnercoords,AllABDPreSynapticLocaions);
[m,n] = ind2sub(size(distMat),find(distMat==0)); % find the abd synapses


PrePostRange_ABD = sqrt(sum((prePartnercoords - postPartnercoords).^2,2))./1000;



% at each of the above motor synapses locations, find inter neurons that
% are inside a certain radius

    prePartnercoords = prePartnercoords.*[5,5,45] ; % convert to nm;
    AllABDiPreSynapticLoactions_nm = AllABDiPostSynapticLocations.*[5,5,45]; % convert to nm 
    ABDidist = [];
for i = 1:length(m) 
    temp = pdist2(prePartnercoords(m(i),:),AllABDiPreSynapticLoactions_nm);
    temp = temp; % covert to um
    ABDidist = [ABDidist;temp];
    %numberOfLocs(i) = length(find(temp<));
    clear temp;
end

%histogram(ABDidist./1000)
%% inters

AllABDiPreSynapticLocaions = [vertcat(ABDIr.PreSynCoords); vertcat(ABDIc.PreSynCoords)];
AllABDPSDlocs = [vertcat(ABDr.PSDID); vertcat(ABDr.PSDID)];
AllABDPostSynapticLocations = PostPartnerCoordinates(AllABDPSDlocs(1:100:end),df);


prePartnercoords = [];
postPartnercoords = [];

for i = 1:size(ABDiPutativeSaccadic.cellIDs,2)
    [partnerID,partnerPSD] = SynapticPartners(ABDiPutativeSaccadic.cellIDs(i),2,df);
    prePartnercoords = [prePartnercoords;PrePartnerCoordinates(partnerPSD,df)];
    postPartnercoords = [postPartnercoords;PostPartnerCoordinates(partnerPSD,df)];
    clear partnerID;
    clear partnerPSD;
end

distMat = pdist2(prePartnercoords,AllABDiPreSynapticLocaions);
[m,n] = ind2sub(size(distMat),find(distMat==0)); % find the abd synapses

PrePostRange_ABDi = sqrt(sum((prePartnercoords - postPartnercoords).^2,2))./1000;



% at each of the above motor synapses locations, find inter neurons that
% are inside a certain radius

    prePartnercoords = prePartnercoords.*[5,5,45] ; % convert to nm;
    AllABDPostSynapticLocations_nm = AllABDPostSynapticLocations.*[5,5,45]; % convert to nm 
    ABDdist = [];
for i = 1:length(m) 
    temp = pdist2(prePartnercoords(m(i),:),AllABDPostSynapticLocations_nm);
    temp = temp; % covert to um
    ABDdist = [ABDdist;temp];
    %numberOfLocs(i) = length(find(temp<));
    clear temp;
end

%%

subplot(1,2,1)
histogram(PrePostRange_ABD)
hold on
histogram(PrePostRange_ABDi)
legend({'EBNm-->ABD','EBNi-->ABDi'});


subplot(1,2,2)

histogram(ABDidist(ABDidist)./1000)
hold on
histogram(ABDdist(ABDdist<5000)./1000)
legend({'EBNm-->ABDi','EBNi-->ABDm'})


