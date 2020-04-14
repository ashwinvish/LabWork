% connectome for Sebastian

clear all;

LoadDataFrame

AllCells = [];
AllCells = unique([unique(df.postsyn_segid(df.postsyn_segid<1e5)); unique(df.presyn_segid(df.presyn_segid<1e5))]);
AllCells = AllCells(AllCells>1e4 & AllCells<1e5);

ConnMatrixPre = zeros(size(AllCells,1));

for i = 1:size(AllCells,1) % pre
        tempPrePartner =  df.presyn_segid(df.postsyn_segid==AllCells(i));
        tempPrePartner = tempPrePartner(tempPrePartner>1e4 & tempPrePartner<1e5);
        if ~tempPrePartner == 0
           ConnMatrixPre(i,:) = histc(tempPrePartner, AllCells);
        end
        ConnMatrixPre(i,i) = 0; % set autapses to 0
        clear tempPrePartner;      
end

 save('AllCells.mat', 'AllCells');
 save('ConnMatrixPre.mat','ConnMatrixPre');

%%

load ConnMatrixPre.mat
load AllCells.mat


load SaccadicProjectingToABDexclusively.mat
load SaccadicProjectingToABDiexclusively.mat
load SaccadicProjectingToABD_ABDi.mat
load ABDPutativeSaccadic.mat
load ABDiPutativeSaccadic.mat
load ABDr.mat
load ABDc.mat
load ABDIr.mat
load ABDIc.mat
load ALX.mat
load DBX.mat

LoadDataFrame;

cellIDs.Saccadic_M = ABDPutativeSaccadic.cellIDs;
cellIDs.Saccadic_I = ABDiPutativeSaccadic.cellIDs;

cellIDs.IBN = [77231,77942,80971,79053,79083,78685,77137,77247,77157,77135,78557,78550,77153,77125,79084,78567,77941,77128];
cellIDs.IBNmirrorPop = [77303,79799,80194,80847,80701,77355,77353,80737,80816,80707,81169,77855,79799,80700,80868];

cellIDs.RS = [79961,81410,77449,77456,77441,77265,77267,77268,78577,77694,77695,78566,80327,82220,78940,82218,...
    82217,76562,77931,77259,77260,79395,79244];


cellIDs.Vestibular_DO = [76631,76700,76656,76655,76660,76679,77393,77395,77375,...
     77815,77803,77807,80591,80579,77255,77697,77364,80223,81126,80760,81117,...
     80672,80622,80572,80569,81122,78426,81328,81575,81582,81870,82161,82165,81553,82169];
 
cellIDs.Vestibular_MO = [76664 76784 76783 77394 77326 77772 77766 77756 77753 77767 77775 ...
     77780 80883 80247 80669 80698 80762 80748 80696 80761 80749 81070 81044 ...
     80991 80967 81008 80883 81045 80986 81023 78647 78619 80247 81141];
 
cellIDs.Vestibular_TO = [78046,78049,78048,78047,78050,78051,78045,80345];



cellIDs.Integrator_r456I = SaccadicProjectingToABDiexclusively.cellID;
cellIDs.Integrator_r456M = SaccadicProjectingToABDexclusively.cellID;
cellIDs.Integrator_r456MI = SaccadicProjectingToABD_ABDi.cellID;
cellIDs.Integrator_r56M = [80163 80167 80177 80179 76688 80204 80206 80210 76682];
cellIDs.Integrator_r78ipis  = temp;
cellIDs.Integrator_r78contra = temp2;

cellIDs.Abducens_M = [[ABDr.cellID],[ABDc.cellID]];
cellIDs.Abducens_I = [[ABDIr.cellID],[ABDIc.cellID]];

% 
cellIDsforMirroring.r78Contra  = [76188,76200,76838,76877,77605,77582,76289,76182,76186,76828,76189,76183,76829,76191,79040,76542,76826,76832,76199,76185,76399];
cellIDsforMirroring.r78ContraPartners  =   [78687,78651,76666,80219,76666,76666,78677,79950,77869,80219,79134,78923,76675,78903,81682,80241,80248,78615,77773,78615,80242];

% save('cellIDs.mat','cellIDs');
% save('cellIDSforMirroring.mat','cellIDsforMirroring');

%% Connectome
A = struct2array(cellIDs);
[~,Matindex] = ismember(A,AllCells);

cspy(ConnMatrixPre(Matindex,Matindex),'Colormap',colorcet('R3','N',15),'Levels',15,'MarkerSize',7);
