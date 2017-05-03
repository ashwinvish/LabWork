function [WD] = UtilityScript_ReadInfoForAllTilesInMontageDir_AV(MontageDirInput)


if nargin < 1
    IsCalledAsSubFunction = false;
else
    IsCalledAsSubFunction = true;
end

if ~IsCalledAsSubFunction
    MontageStackDir = uigetdir('Z:\Hayworth\MasterUTSLDir\Cortex\w007\MontageStack_08\','Pick the montage directory');
    if MontageStackDir == 0
        return;
    end
else
    MontageStackDir = MontageDirInput;
end

MaxR = 1;
MaxC = 1;

n=1;

DirAndWildCardStr = sprintf('%s/*_Montage',MontageStackDir);
ListOfMontageDirectories = dir(DirAndWildCardStr);

OrderedArrayOfMontageDirectoryNames = [];

for i = 1:length(ListOfMontageDirectories)
    TempArray = strfind(ListOfMontageDirectories(i).name,'_');
    SectionNumStr = ListOfMontageDirectories(i).name(TempArray(end-1)+4:TempArray(end)-1);
    SectionNum = str2num(SectionNumStr);
    
    OrderedCellArrayOfMontageDirectoryNames{SectionNum} = ListOfMontageDirectories(i).name;
    
end

MaxR = 1;
MaxC = 1;

ArrayOfTileInfos = repmat(struct('FileNameAndPath',[],'R' , [] , 'C' , [], 'Info' ,[]),1,length(OrderedCellArrayOfMontageDirectoryNames));


MovieFrameNum = 0;

for SectionNum = 1:length(OrderedCellArrayOfMontageDirectoryNames)
    MontageDirName = OrderedCellArrayOfMontageDirectoryNames{SectionNum};
    
    if ~isempty(MontageDirName)
        MovieFrameNum = MovieFrameNum + 1;
        MontageDir = sprintf('%s/%s', MontageStackDir, MontageDirName);
        
        ListOfTifFiles = dir(MontageDir);
        
        for i = 1:length(ListOfTifFiles)
            FileName = ListOfTifFiles(i).name;
            if length(FileName) >= 5
                if strcmp('Tile_',FileName(1:5))
                    %disp(sprintf('Loading file: %s',FileName));
                    FileNameAndPath = sprintf('%s/%s', MontageDir, FileName);
                    
                    
                    
                    %Tile_r1-c6_AW01_sec65
                    if strcmp('.mat', FileName((end-3):end))
                        TempArray = strfind(FileName,'_');
                        SubStr = FileName(TempArray(1)+2:end);
                        TempArray2 = strfind(SubStr,'-');
                        RowStr = SubStr(1:TempArray2(1)-1);
                        TempArray3 = strfind(SubStr,'_');
                        ColStr = SubStr(TempArray2(1)+2:TempArray3(1)-1);
                        
                        
                        R = str2num(RowStr);
                        C = str2num(ColStr);
                        
                        if MaxR < R
                            MaxR = R;
                        end
                        if MaxC < C
                            MaxC = C;
                        end
                        
                        if ~IsCalledAsSubFunction
                            disp(sprintf('File: %s, (r,c) = (%d, %d)',FileName, R, C));
                        end
                        
                        ArrayOfTileInfos(n).FileNameAndPath = FileNameAndPath;
                        ArrayOfTileInfos(n).R = R;
                        ArrayOfTileInfos(n).C = C;
                        
                        
                        
                        n = n+1;
                    end
                end
            end
        end
    end
end



WD = zeros(length(ArrayOfTileInfos),1);
for n = 1:length(ArrayOfTileInfos)
    %INFO = IMFINFO(FILENAME,FMT)
    
    
    R = ArrayOfTileInfos(n).R;
    C = ArrayOfTileInfos(n).C;
    
   % disp(sprintf(' Loading File: %s, (r,c) = (%d, %d)',ArrayOfTileInfos(n).FileNameAndPath, R, C));
    load(ArrayOfTileInfos(n).FileNameAndPath, 'Info');
    ArrayOfTileInfos(n).Info = Info;
    WD(n) = Info.WorkingDistance*1000 ; 
    %disp(sprintf(' Info.WorkingDistance = %0.5g', Info.WorkingDistance));
   
end


