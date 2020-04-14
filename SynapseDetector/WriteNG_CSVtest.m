function WriteNG_CSVtest(pointCloud,cellID)

filepath = '/Users/ashwin/Google Drive/Zfish/IOcoordinates';
filename = sprintf('%s/%d_%s.csv',filepath,cellID,'toRemove');

fileID = fopen(filename,'w');

for i = 1:size(pointCloud,1)
    fprintf(fileID,'%c%c%d%c%d%c%d%c%c%c%c%c%c%c%c%c%s%c\n',...
        '"','(',pointCloud(i,1),',',pointCloud(i,2),',',pointCloud(i,3),')','"',',',',',',',',',',',',',',','Point',',');
end
fclose(fileID);

end

