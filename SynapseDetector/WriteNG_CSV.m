function [inputSize,outputSize] = WriteNG_CSV(cellID,df)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[~,inputs] = SynapticPartners(cellID,1,df); % on the dendrite
[~,outputs] = SynapticPartners(cellID,2,df); % on the axon

inputSize = length(inputs);
outputSize = length(outputs);

filepath = '/Users/ashwin/Google Drive/Zfish/IOcoordinates';

input_coords = PrePartnerCoordinates(inputs,df);
output_coords = PostPartnerCoordinates(outputs,df);

if inputSize>0
    input_xcoord = input_coords(:,1);
    input_ycoord = input_coords(:,2);
    input_zcoord = input_coords(:,3);
    
    % input CSV file
    filename = sprintf('%s/%d_%s.csv',filepath,cellID,'inputs');
    fileID = fopen(filename,'w');
    for i = 1:size(input_xcoord,1)
        fprintf(fileID,'%c%c%d%c%d%c%d%c%c%c%c%c%c%c%c%c%s%c\n',...
            '"','(',input_xcoord(i),',',input_ycoord(i),',',input_zcoord(i),')','"',',',',',',',',',',',',',',','Point',',');
    end
    fclose(fileID);
    
end

if outputSize>0
    output_xcoord = output_coords(:,1);
    output_ycoord = output_coords(:,2);
    output_zcoord = output_coords(:,3);
    
    % output CSV file
    filename = sprintf('%s/%d_%s.csv',filepath,cellID,'outputs');
    fileID = fopen(filename,'w');
    for i = 1:size(output_xcoord,1)
        fprintf(fileID,'%c%c%d%c%d%c%d%c%c%c%c%c%c%c%c%c%s%c\n',...
            '"','(',output_xcoord(i),',',output_ycoord(i),',',output_zcoord(i),')','"',',',',',',',',',',',',',',','Point',',');
    end
    fclose(fileID);
end

filename = fullfile(filepath,'io.txt');
fileID = fopen(filename,'a')';
fprintf(fileID,' %d %d %d\n', cellID, inputSize, outputSize);
fclose(fileID);

end

