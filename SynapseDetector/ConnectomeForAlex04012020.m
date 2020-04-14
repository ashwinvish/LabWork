% connectome for Alex - 03302020

clc;
clear;

load block1_subMod.mat % onto ABDi
load block2_subMod.mat % onto ABD

LoadDataFrame;

load AllCells.mat
load ConnMatrixPre.mat

ABDr_CellIDs = [82140,82145,82143,77648,82146,77302,77301,82194,77705,77710,77305,77672,77300,77709,77661];
ABDc_CellIDs = [82213,77654,77646,82212,77295,81172,77657,82197,77682,77154, 77652,77658,77628,77292,77688,82195,77296];
ABDIr_CellIDs = [78553,77150,77886,78547,77631,77158,77665,78552,77668,77618,77634,78556];
ABDIc_CellIDs =[78574,79051,77148,77625,77144,77692,77641,79066,77640,77643];

cellIDType_04012020 = [repmat('Int2a',length(block1_subMod_gamm075.cellID),1);...
              repmat('Int2b',length(block2_subMod_gamm075.cellID),1);...
              repmat('ABD_m',length(ABDr_CellIDs)+length(ABDc_CellIDs),1);...
              repmat('ABD_i',length(ABDIr_CellIDs)+length(ABDIc_CellIDs),1)];

A = [block1_subMod_gamm075.cellID;...
           block2_subMod_gamm075.cellID;...
           ABDr_CellIDs';ABDc_CellIDs';ABDIr_CellIDs';ABDIc_CellIDs'];
       
for i = 1:size(A,1)
    tempA = SynapticPartners(A(i),1,df);
    tempB = tempA(tempA~=A(i)); % remove autpses
    totalInputs_04012020(i) = length(tempB);
    clear tempA;
    clear tempB;
end
% 
[~,Matindex] = ismember(A,AllCells);

connMatrixforAlex_04012020 = ConnMatrixPre(Matindex,Matindex);
cspy(connMatrixforAlex_04012020,'Colormap',colorcet('R3','N',15),'Levels',15,'MarkerSize',7);

% save('cellIDType_04012020.mat','cellIDType_04012020');
% save('totalInputs_04012020.mat','totalInputs_04012020');
% save('connMatrixforAlex_04012020','connMatrixforAlex_04012020');
