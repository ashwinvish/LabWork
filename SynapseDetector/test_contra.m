% read contraCell list
contraList = dlmread('/usr/people/ashwinv/seungmount/research/Ashwin/SynapseDetector/contraList-07032018.txt','\t',1,0);
contra = [];
ipsi = [];
for i = 1:size(contraList,1)
    if contraList(i,3) == 1
    contra(i) = contraList(i,1);
    else 
     ipsi(i) = contraList(i,1);
    end
end

functionalCellIDs_old = [76198, 76199, 76200, 76201, 76181, 76187, 76210, 76182, 76183, 76185, 76186, 76189, 76184, 76190, ...
    76191, 76192, 76193, 76194, 76195, 76196,76188, 76197];
functionalCellIDs_new = [76181,76182,76183,76184,76185,76186,76187,76188,76189,76190,76191 ...
    ,76192,76193,76194,76195,76196,76197,76198,76199,76200,76201,76210];

% old and new cell sort order cellIds_new = cellIds_old(sortOrder)
sortOrder = [];
for i = 1:1:size(functionalCellIDs_new,2)
    sortOrder = [sortOrder, find(functionalCellIDs_old==functionalCellIDs_new(i))];
end

% classify presynaptic partners into ipsi and contra
for ii = 1:size(functionalCellIDs_new,2)
    A = [];
    B = [];
    m = [];
    [A,B] = SynapticPartners(functionalCellIDs_new(ii),1,df);
    m = find(A<1e5);
    preSynapsePsdID = B(m);
    preSynapseCellID = A(m);
    
    contraPreSynapseCoord = [];
    ipsiPreSynapseCoord = [];
    contraPreSynapseCellID = [];
    ipsiPreSynapseCellID = [];

    for i = 1:size(preSynapseCellID,1)
        if ismember(preSynapseCellID(i),contra)
            contraPreSynapseCoord = [contraPreSynapseCoord; PrePartnerCoordinates(preSynapsePsdID(i),df)];
            contraPreSynapseCellID = [contraPreSynapseCellID;preSynapseCellID(i)];
        else
            ipsiPreSynapseCoord = [ipsiPreSynapseCoord; PrePartnerCoordinates(preSynapsePsdID(i),df)];
            ipsiPreSynapseCellID = [ipsiPreSynapseCellID, preSynapseCellID(i)];
        end
    end
    
    contraIntPreSynapseCellID{ii} = contraPreSynapseCellID;
    ipsiIntPreSynapseCellID{ii} = ipsiPreSynapseCellID;
    numberOfContraPreSynapses(ii) = size(contraPreSynapseCoord,1);
    numberOfIpsiPreSynapses(ii) = size(ipsiPreSynapseCoord,1);
end

% figures
CellDisplayOrder = [21,1,7,4,12,17,19,20,2,3,5,6,9,11,8,18,10,13,14,15,16,22];
figure; 
subplot (4,7,6);
plot(1:22 , numberOfIpsiPreSynapses(CellDisplayOrder)./numberOfContraPreSynapses(CellDisplayOrder),'.','MarkerSize',18);
hold on;
line([6.5,6.5],[0,10],'Color',[0.5,0.5,0.5],'LineStyle','--');
line([15.5,15.5],[0,10],'Color',[0.5,0.5,0.5],'LineStyle','--');
line([21.5,21.5],[0,10],'Color',[0.5,0.5,0.5],'LineStyle','--');
set(gca,'XLim',[0,23]);
xlabel('Cell number')
ylabel({'Ipsi. to Contra inputs';'~ E to I ratio'});
%axis square;
box off;

%%
AbdCells = [77646,77296,77628,77652,77292,77688,77682,77658,77654,77662,77657,77295,77648,77710,77705,77305,77301,77672,77709,77661,77300,77302,77641,77643,77144,77625,77692,77640,79066,79051,77148,77668,77631,78552,77665,77618,77150,77886,78547,77158,78556,77634];
contraProjectingToInt = [];
ipsiProjectingToInt = [];
for i = 1:22
    contraProjectingToInt = [contraProjectingToInt;contraIntPreSynapseCellID{i}];
    ipsiProjectingToInt = [ipsiProjectingToInt, ipsiIntPreSynapseCellID{i}];
end

contraProjectingToAbd = [];
ipsiProjectingToAbd = [];

for i = 1:size(AbdCells,2)
    A = [];
    B = [];
    [A,B] = SynapticPartners(AbdCells(i),1,df);
    m = find(A<1e5);
    B = B(m); % presynapse PSD ID
    A = A(m); % presynapse cellID
    C = intersect(A,contraProjectingToInt);
    D = intersect(A, ipsiProjectingToInt);
    contraProjectingToAbd = [contraProjectingToAbd; C];
    ipsiProjectingToAbd = [ipsiProjectingToAbd; D];
end