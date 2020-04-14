% find EBNs
clear;
LoadDataFrame

load ABDr.mat
load ABDc.mat
load ABDIr.mat
load ABDIc.mat

load IntBlock.mat
load TurnBlock.mat

IBNall = [82714,82710,79083,82264,78550,79084,82712,82709,82713,83184,77125,77942,77231,77247,77153,78685,77137,83183];

allABDm.Inputs = [vertcat(ABDr.Inputs);vertcat(ABDc.Inputs)];
allABDm.Inputs = unique(allABDm.Inputs);
allABDm.Inputs = allABDm.Inputs(allABDm.Inputs<1e5);

allABDi.Inputs = [vertcat(ABDIr.Inputs);vertcat(ABDIc.Inputs)];
allABDi.Inputs = unique(allABDi.Inputs);
allABDi.Inputs = allABDi.Inputs(allABDi.Inputs<1e5);

EBNm = [];
for i = 1:length(allABDm.Inputs)
    A = SynapticPartners(allABDm.Inputs(i),2,df);
    if (sum(ismember(A,IBNall)) & sum(ismember(A,Block3)) & sum(isPostSynapseIntegrator(A,df)))>0
        EBNm = [EBNm;allABDm.Inputs(i)];
    end
    clear A;
end

EBNi = [];
for i = 1:length(allABDi.Inputs)
    A = SynapticPartners(allABDi.Inputs(i),2,df);
    if (sum(ismember(A,IBNall)) & sum(ismember(A,Block3)) & sum(isPostSynapseIntegrator(A,df)))>0
        EBNi = [EBNi;allABDi.Inputs(i)];
    end
    clear A;
end

EBNm  = setdiff(EBNm,Block3);
EBNi = setdiff(EBNi,Block3);

% common partners

EBNm_EBNi = intersect(EBNm,EBNi);
EBNmonly = setdiff(EBNm,EBNm_EBNi);
EBNionly = setdiff(EBNi,EBNm_EBNi);

EBNmonly_cleaned = [80829,78836,81839,79712,81868,79728,79838,79724,80841,79872];
