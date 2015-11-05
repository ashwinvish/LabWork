clc;
clear all;

folderName = uigetdir;
%d=dir(folder_name);
NumSections = 164;
row = 2; % number of rows in the reimaged folder
column = 4; % number of columns in the reimaged folder
sec = 1; % starting section number
wafernumber = 9; % wafer number
targetTemplateTif = 'Tile_r%d-c%d_W%03d_sec%d.tif';
targetTemplateMat = 'Tile_r%d-c%d_W%03d_sec%d.mat';



for kk = sec:1:NumSections

    for i = 1:1:row
        for j = 1:1:column

            newRow = i+5; % i+6 , new row number
            newCol = j+2; % j+2 , new column number

            if (i==1 && j==4); continue; end;

            secDir = sprintf('W%03d_Sec%d_Montage',wafernumber,kk);

            OldFnameTif = sprintf(targetTemplateTif,i,j,wafernumber,kk);
            OldFnameMat = sprintf(targetTemplateMat,i,j,wafernumber,kk);

            path1Tif = fullfile(folderName,secDir,OldFnameTif);
            path1Mat = fullfile(folderName,secDir,OldFnameMat);

            NewFnameTif = sprintf(targetTemplateTif,newRow,newCol,wafernumber,kk);
            NewFnameMat = sprintf(targetTemplateMat,newRow,newCol,wafernumber,kk);

            path2Tif = fullfile(folderName,secDir,NewFnameTif);
            path2Mat = fullfile(folderName,secDir,NewFnameMat);

            movefile (path1Tif , path2Tif);
            disp(sprintf('Moved %s to %s' ,OldFnameTif ,NewFnameTif ));
            movefile (path1Mat , path2Mat);
            disp(sprintf('Moved %s to %s' ,OldFnameMat ,NewFnameMat ));
            %delete(path1Tif, path1Mat);

        end
    end
end
