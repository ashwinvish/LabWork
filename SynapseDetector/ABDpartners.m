clc;

if ismac
    addpath(genpath('/Users/admin/Documents/Scripts'));
    df = readtable('/Users/admin/Documents/SynapseDetector/140318_134755_final.csv');
else
    addpath(genpath('/usr/people/ashwinv/seungmount/research/Ashwin/Scripts'));
    df = readtable('/usr/people/ashwinv/seungmount/research/Ashwin/SynapseDetector/09202018.csv');
end

ABDr_CellIDs = [77648, 77710, 77300, 77705, 77305, 77301, 77709, 77672, 77302];
ABDc_CellIDs = [77154, 77646, 77682 ,77628 ,77295 , 77652 ,77292 ,77688 ,77654 ,77658 ,77657 ,77662, 77296];
ABDIr_CellIDs = [77631, 77150, 77618, 77886, 78547, 77158, 78556, 78552, 77665, 77668, 77634, 78553];
ABDIc_CellIDs = [77148, 77625, 77641, 77692, 77144, 77643, 77640, 79051, 79066, 78574];

ABD_All = [ABDr_CellIDs,ABDc_CellIDs,ABDIr_CellIDs,ABDIc_CellIDs];

ABDr_inputs = inputsOfClass(ABDr_CellIDs,df);
ABDr_inputs(ABDr_inputs==0) = [];
ABDr_inputs = sort(ABDr_inputs);

ABDc_inputs = inputsOfClass(ABDc_CellIDs,df);
ABDc_inputs(ABDc_inputs==0) = [];
ABDc_inputs = sort(ABDc_inputs);

ABDIr_inputs = inputsOfClass(ABDIr_CellIDs,df);
ABDIr_inputs(ABDIr_inputs==0) = [];
ABDIr_inputs = sort(ABDIr_inputs);

ABDIc_inputs = inputsOfClass(ABDIc_CellIDs,df);
ABDIc_inputs(ABDIc_inputs==0)= [];
ABDIc_inputs = sort(ABDIc_inputs);


ABDcommon = intersect(ABDr_inputs,ABDc_inputs);
ABDIcommon = intersect(ABDIr_inputs,ABDIc_inputs);

ABDrOnly = setdiff(ABDr_inputs,ABDcommon);
ABDcOnly = setdiff(ABDc_inputs,ABDcommon);

ABDIrOnly = setdiff(ABDIr_inputs,ABDIcommon);
ABDIcOnly = setdiff(ABDIc_inputs,ABDIcommon);



addpath(genpath('../')); 
ml = readtable('/usr/people/ashwinv/seungmount/research/Ashwin/SynapseDetector/ManualSortList-102218.csv');
mlMat = table2array(ml);
load('IntConnMatrix.mat')
load('IntPartners.mat')
load('AllCells.mat');
load('ConnMatrixPre.mat');

for i = 1:size(ABD_All,2)
   tempID = find(ABD_All(i) == AllCells);
   ABDmatrix(i,:) = ConnMatrixPre(tempID,:);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





