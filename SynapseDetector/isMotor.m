function [] = isMotor(cellID)

ABDr_CellIDs = [77648, 77710, 77300, 77705, 77305, 77301, 77709, 77672, 77302];
ABDc_CellIDs = [77154, 77646, 77682 ,77628 ,77295 , 77652 ,77292 ,77688 ,77654 ,77658 ,77657 ,77662, 77296];
ABDIr_CellIDs = [77631, 77150, 77618, 77886, 78547, 77158, 78556, 78552, 77665, 77668, 77634, 78553];
ABDIc_CellIDs = [77148, 77625, 77641, 77692, 77144, 77643, 77640, 79051, 79066, 78574];

if ismember(cellID,ABDr_CellIDs)
    logic = 1;
    type = char('ABDr');
elseif ismember (cellID, ABDc_CellIDs)
    logic = 1;
    type = char('ABDc');
elseif ismember (cellID, ABDIr_CellIDs)
    logic = 1;
    type = char('ABDIr');
elseif ismember(cellID,ABDIc_CellIDs)
    logic = 1;
    type = char('ABDIc');
else
    logic = 0
end

end