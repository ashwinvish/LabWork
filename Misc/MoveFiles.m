% script to move files from one folder to another
% created by Ashwin Vishwanathan-03/20/2015
% contact ashwinv@princeton.edu

clc;
clear all;

% get folder

%originalFolder = uigetdir();
originalFolder = ('/local_data/Talmo/Zfish/W009/HighResImages_fine_5nm_120apa_W009_2ndPass/');
%NewFolder = uigetdir();
NewFolder = ('/local_data/Talmo/Zfish/W009/');

OriginalFolderAttributes = dir(originalFolder);


for i = 1:length(OriginalFolderAttributes)
    cd (originalFolder);
    if isempty(strfind(OriginalFolderAttributes(i).name,'Montage'))==1
        continue
    else
        cd (OriginalFolderAttributes(i).name);
        movefile('Tile*.tif', fullfile(NewFolder,OriginalFolderAttributes(i).name));
        str = sprintf('Finished moving files in %s',OriginalFolderAttributes(i).name);
        disp(str);
    end
end



        
        

