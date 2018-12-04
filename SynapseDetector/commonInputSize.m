function[commSize1, commSize2] = commonInputSize(cellID1,cellID2, commonCellID, df)

commSize1 = [];
for i=1:size(commonCellID,1)
    [A,B] = SynapticPartners(cellID1,1,df);
    commonInputLocation = find(A==commonCellID(i));
    commonInputPSDID =B(commonInputLocation);
    commSize1 = [commSize1;df.size(commonInputPSDID)];
end
commSize1= sum(commSize1);

commSize2 = [];
for i=1:size(commonCellID,1)
    [A,B] = SynapticPartners(cellID2,1,df);
    commonInputLocation = find(A==commonCellID(i));
    commonInputPSDID =B(commonInputLocation);
    commSize2 = [commSize2;df.size(commonInputPSDID)];
end
commSize2 = sum(commSize2);
end



