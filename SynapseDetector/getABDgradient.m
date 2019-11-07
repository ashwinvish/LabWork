function [TargetQueryPathlength,TargetGradient,TargetingABDTotalInputs] = getABDgradient(TargetingABDCellIDs,QueringCellIDs,isNormalized,withAIS)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

if nargin<4
    withAIS = false
end

if nargin <3
    isNormalized = false;
end

load MotorAIS.mat

if ismac
    addpath(genpath('/Users/ashwin/Documents/'));
    df = readtable('/Users/ashwin/Google Drive/Zfish/SynapseDetector/04152019.csv');
else
    addpath(genpath('/usr/people/ashwinv/seungmount/research/Ashwin/Scripts'));
    df = readtable('/usr/people/ashwinv/seungmount/research/Ashwin/SynapseDetector/11252018.csv');
end



if strcmp(TargetingABDCellIDs,'ABDr')
    load ABDr.mat
    TargetingABDCellIDs = ABDr;
    for i = 1:numel(ABDr)
        TargetingABDCellIDs(i).motorAIStransformed = TransformPoints(MotorAIS(i,:),0);
    end
    
else if strcmp(TargetingABDCellIDs,'ABDc')
        load ABDc.mat
        TargetingABDCellIDs = ABDc;
        for i = 1:numel(ABDc)
            TargetingABDCellIDs(i).motorAIStransformed = TransformPoints(MotorAIS(14+i,:),0);
        end
    end
end

%TargetQueryPathlength = cell(numel(TargetingABDCellIDs));

for i = 1:numel(TargetingABDCellIDs)
    if ~isempty(TargetingABDCellIDs(i).Tree)
        if withAIS
            %TargetingABDTotalInputs(i) = numel(TargetingABDCellIDs(i).Inputs);
            temp1 = ismember(TargetingABDCellIDs(i).Inputs,QueringCellIDs);
            TargetQueryPathlength{i} = TargetingABDCellIDs(i).PathLength(temp1)'./(max(Pvec_tree(TargetingABDCellIDs(i).Tree{1}))-PathLengthToCoordinate(TargetingABDCellIDs(i).motorAIStransformed,TargetingABDCellIDs(i).Tree{1}));
        else
            temp1 = ismember(TargetingABDCellIDs(i).Inputs,QueringCellIDs);
            TargetQueryPathlength{i} = TargetingABDCellIDs(i).PathLength(temp1)'./max(Pvec_tree(TargetingABDCellIDs(i).Tree{1}));
        end
    end
    clear temp1;
end

%TargetGradient = zeros(numel(TargetingABDCellIDs),10);

for i = 1:numel(TargetingABDCellIDs)
    if  ~isempty(TargetingABDCellIDs(i).Tree)
        TargetGradient(i,:) =  histcounts(TargetQueryPathlength{i},0:0.1:1);
        if isNormalized
            TargetGradient(i,:) =  histcounts(TargetQueryPathlength{i},0:0.1:1,'Normalization','probability');
        end
    end
end


end

