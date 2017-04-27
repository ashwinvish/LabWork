function [ SigNoise ] = SNRall( input_args )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% read part of image
% if ~IsCalledAsSubFunction
    MontageStackDir = uigetdir('/gaba/atlas/MasterUTSLdirectory/','Pick the montage stack directory');
    if MontageStackDir == 0
        return;
    end
% else
%     MontageDir = MontageDirInput;
% end

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

MovieFrameNum = 0;n=1;
for SectionNum = 1:length(OrderedCellArrayOfMontageDirectoryNames)
    MontageDirName = OrderedCellArrayOfMontageDirectoryNames{SectionNum};
    
    if ~isempty(MontageDirName)
        MovieFrameNum = MovieFrameNum + 1;
        MontageDir = sprintf('%s/%s', MontageStackDir, MontageDirName);
        % disp(sprintf('Analizing: %s',  MontageDir));
        
        ListOfTifFiles = dir(MontageDir);
        for i = 1:length(ListOfTifFiles)
            FileName = ListOfTifFiles(i).name;
            if length(FileName) >= 5
                if strcmp('Tile_',FileName(1:5))
                    %disp(sprintf('Loading file: %s',FileName));
                    FileNameAndPath = sprintf('%s/%s', MontageDir, FileName);
                    
                    %Tile_r1-c6_AW01_sec65
                    if ~strcmp('.mat', FileName((end-3):end))
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

%
clear i;
SigNoise = zeros(length(ArrayOfTileInfos),1);

for i = 1: length(ArrayOfTileInfos)
    fname = ArrayOfTileInfos(i).FileNameAndPath;
    im1 = imread(fname);
    disp(sprintf('Loading File: %s',fname));
    SigNoise(i) = SNR(im1);
    i
end
    
    %
    
 k =0;
 clear i;
 
 % plot SNR for all sections
 
    for i=1:169
        plot(i,SigNoise(1+k:16+k),'o');
        hold all;
        k = k+16;
    end
    
    xlabel('Section Number');
    ylabel('SNR');
    
    figure;
    
    hist(SigNoise(:,1),100);
    xlabel('SNR');
    ylabel('Frequency');
    mu = mean(SigNoise(:,1));
    hold on;
    plot([mu,mu],ylim, 'r--','LineWidth' ,2);
    hold off;
    
end


