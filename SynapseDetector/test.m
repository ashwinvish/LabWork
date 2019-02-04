for i = 1:numel(ABDr_CellIDs)
    ABDr(i) = InputsByClass(ABDr_CellIDs(i),df);
    filename = sprintf('%d.swc',ABDr_CellIDs(i));
    if exist(fullfile(fname,filename))>0
        temp = dlmread(fullfile(fname,filename), ' ');
        coord = temp(:,3:5);
        ABDr(i).Origin = TransformPoints(coord(1,:),4); 
    else 
        ABRr(i).Origin = NaN;
    end
    clear coord;
    clear temp;
end