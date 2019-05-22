function [yesNoCell] = isPostSynapseMotor(cellID,df)
%   isPostSynapseMotor returns a Nx4 vector containing the number of synapses
% onto each of the abducens nuclei
%   cellID is an vector with the list of cells ID that are being queried.

if isempty(cellID)
    yesNoCell = [];
end

ABDr_CellIDs = [82140,82145,82143,77648,82146,77302,77301, 82194,77705,77710,77305,77672,77300,77709];
ABDc_CellIDs = [82213,77654,77646,82212,77295,81172,77657,82197,77682,77154, 77652,77658,77628,77292,77688,82195,77296];
ABDIr_CellIDs = [78553,77150,77886,78547,77631,77158,77665,78552,77668,77618,77634];
ABDIc_CellIDs = [78574,79051,77148,77625,77144,77692,77641,79066,77640,77643];


for jj = 1:numel(cellID)
    [A,~] = SynapticPartners(cellID(jj),2,df); % post synaptic partners of axon
    A = A(A<1e5);
   % A = unique(A);
    
    for i = 1:length(A)
        yesNo(1,1:numel(ABDr_CellIDs),i) = ismember(ABDr_CellIDs,A(i));
        yesNo(2,1:numel(ABDc_CellIDs),i) = ismember(ABDc_CellIDs,A(i));
        yesNo(3,1:numel(ABDIr_CellIDs),i) = ismember(ABDIr_CellIDs,A(i));
        yesNo(4,1:numel(ABDIc_CellIDs),i) = ismember(ABDIc_CellIDs,A(i));
    end
    
    %yesNoCell(jj,:) = sum(yesNo,2);
    yesNoCell = sum(yesNo(:,:,:),3);
    
end

end


