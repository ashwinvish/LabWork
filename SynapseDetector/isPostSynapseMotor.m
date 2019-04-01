function [yesNoCell] = isPostSynapseMotor(cellID,df)
%   isPostSynapseMotor returns a Nx4 vector containing the number of synapses
% onto each of the abducens nuclei
%   cellID is an vector with the list of cells ID that are being queried.

if isempty(cellID)
    yesNoCell = [];
end

ABDr_CellIDs = [77648, 77710, 77300, 77705, 77305, 77301, 77709, 77672, 77302];
ABDc_CellIDs = [77154, 77646, 77682 ,77628 ,77295 , 77652 ,77292 ,77688 ,77654 ,77658 ,77657 ,77662, 77296];
ABDIr_CellIDs = [77631, 77150, 77618, 77886, 78547, 77158, 78556, 78552, 77665, 77668, 77634, 78553];
ABDIc_CellIDs = [77148, 77625, 77641, 77692, 77144, 77643, 77640, 79051, 79066, 78574];


for jj = 1:numel(cellID)
    [A,~] = SynapticPartners(cellID(jj),2,df); % post synaptic partners of axon
    A = A(A<1e5);
    A = unique(A);
    
    for i = 1:length(A)
        yesNo(i,1) = sum(ismember(ABDr_CellIDs,A(i)));
        yesNo(i,2) = sum(ismember(ABDc_CellIDs,A(i)));
        yesNo(i,3) = sum(ismember(ABDIr_CellIDs,A(i)));
        yesNo(i,4) = sum(ismember(ABDIc_CellIDs,A(i)));
    end
    
    yesNoCell(jj,:) = sum(yesNo,1);
end

end


