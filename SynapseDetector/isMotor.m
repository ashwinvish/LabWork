function [edges] = isMotor(cellID,df)
if ~isempty(cellID)
    ABDr_CellIDs = [77648, 77710, 77300, 77705, 77305, 77301, 77709, 77672, 77302];
    ABDc_CellIDs = [77154, 77646, 77682 ,77628 ,77295 , 77652 ,77292 ,77688 ,77654 ,77658 ,77657 ,77662, 77296];
    ABDIr_CellIDs = [77631, 77150, 77618, 77886, 78547, 77158, 78556, 78552, 77665, 77668, 77634, 78553];
    ABDIc_CellIDs = [77148, 77625, 77641, 77692, 77144, 77643, 77640, 79051, 79066, 78574];
    
    ABD_All = [ABDr_CellIDs,ABDc_CellIDs,ABDIr_CellIDs,ABDIc_CellIDs];
    % give each type an ID
    ABD_All_type = [ones(1,size(ABDr_CellIDs,2)),2*ones(1,size(ABDc_CellIDs,2)),...
        3*ones(1,size(ABDIr_CellIDs,2)),4*ones(1,size(ABDIc_CellIDs,2))] ;
    
    %UniqueCells = unique(cellID);
    
    for i =1:size(cellID,1)
        [A,B] = SynapticPartners(cellID(i),2,df);
        [lia, lob] =  ismember(A,ABD_All);
        temp = ABD_All_type(lob(lob>0))';
        edges(i,:) = [cellID(i),histcounts(temp,[0.9,1.9,2.9,3.1,4])];
    end
    
else
    edges = [];
end
end