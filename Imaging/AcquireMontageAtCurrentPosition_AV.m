function AcquireMontageAtCurrentPosition_AV(WaferName, LabelStr)
%This function uses the current scan rotation to determine how to move the
%stage from the current stage position. It assumes that the calling
%function has properly setup the scan rotation.

WaferName = 'S2-W001';
LabelStr = 't76';
global GuiGlobalsStruct ;

LogFile_WriteLine(['Beginning section ' LabelStr])

if exist([GuiGlobalsStruct.TempImagesDirectory '\watchQ.mat'])
    load([GuiGlobalsStruct.TempImagesDirectory '\watchQ.mat'], 'q')
else
    q.fileNum = 0;
end

tilesTaken = {};

disp('Turning stage backlash ON in X and Y');
GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeString('DP_X_BACKLASH','+ -');
GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeString('DP_Y_BACKLASH','+ -');
GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeString('DP_STAGE_BACKLASH', 'On');
pause(.2);

%Get current stage position (this will be center of montage)
StageX_Meters_CenterOfMontage = GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeSingle('AP_STAGE_AT_X');
StageY_Meters_CenterOfMontage = GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeSingle('AP_STAGE_AT_Y');
stage_z = GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeSingle('AP_STAGE_AT_Z');
stage_t = GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeSingle('AP_STAGE_AT_T');
stage_r = GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeSingle('AP_STAGE_AT_R');
stage_m = GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeSingle('AP_STAGE_AT_M');
MyStr = sprintf('(In AcquireMontageAtCurrentPosition) Stage Position(x,y,z,t,r,m) = (%0.7g, %0.7g, %0.7g, %0.7g, %0.7g, %0.7g, )'...
    ,StageX_Meters_CenterOfMontage,StageY_Meters_CenterOfMontage,stage_z,stage_t,stage_r, stage_m);
disp(MyStr);
disp(' ');


RowDistanceBetweenTileCentersInMicrons = GuiGlobalsStruct.MontageTarget.MontageTileWidthInMicrons * ...
    (1-GuiGlobalsStruct.MontageTarget.PercentTileOverlap/100);
ColDistanceBetweenTileCentersInMicrons = GuiGlobalsStruct.MontageTarget.MontageTileHeightInMicrons * ...
    (1-GuiGlobalsStruct.MontageTarget.PercentTileOverlap/100);
NumRowTiles = GuiGlobalsStruct.MontageTarget.NumberOfTileRows;
NumColTiles = GuiGlobalsStruct.MontageTarget.NumberOfTileCols;

%Setup unit vectors in the
theta_Degrees = GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeSingle('AP_SCANROTATION');
theta_Radians = (pi/180)*theta_Degrees;
cosTheta = cos(theta_Radians);
sinTheta = sin(theta_Radians);

c_target_north_UnitVector = sinTheta;
r_target_north_UnitVector = -cosTheta;

c_target_east_UnitVector = cosTheta;
r_target_east_UnitVector = sinTheta;

if (NumRowTiles == 1) && (NumColTiles == 1) %if single image montage then put in main directory
    MontageDirName = GuiGlobalsStruct.TempImagesDirectory;
else
    MontageDirName = sprintf('%s\\%s_Sec%s_Montage', GuiGlobalsStruct.TempImagesDirectory,WaferName, LabelStr);
    if ~exist(MontageDirName,'dir')
        disp(sprintf('Creating directory: %s',MontageDirName));
        [success,message,messageid] = mkdir(MontageDirName);
    end
end

StitchFigNum = 1234;
figure(StitchFigNum);
clf;

%pre allocate StageStitchedImage
ImageWidthInPixels = GuiGlobalsStruct.MontageParameters.TileWidth_pixels;
MaxSSTileR = 256;
Increment = floor(ImageWidthInPixels/MaxSSTileR); %produce 256x256 images
MaxSSTileC = MaxSSTileR; %tiles are square
StageStitchedImage = uint8(255*ones(MaxSSTileR*NumRowTiles, MaxSSTileC*NumColTiles));
DummyTile = 255*ones(MaxSSTileR, MaxSSTileC);
BorderPixels = 3;
DummyTile(1:BorderPixels,:) = 0;
DummyTile(end-BorderPixels+1:end,:) = 0;
DummyTile(:,1:BorderPixels) = 0;
DummyTile(:,end-BorderPixels+1:end) = 0;
for RowIndex = 1:NumRowTiles
    for ColIndex = 1:NumColTiles
        StartR = (MaxSSTileR*(RowIndex-1))+1;
        StartC = (MaxSSTileC*(ColIndex-1))+1;
        StageStitchedImage(StartR:StartR+MaxSSTileR-1, StartC:StartC+MaxSSTileC-1) = DummyTile;
        
        StageStitched_TextStringsArray(RowIndex, ColIndex).textX = 0;
        StageStitched_TextStringsArray(RowIndex, ColIndex).textY = 0;
        StageStitched_TextStringsArray(RowIndex, ColIndex).Text= '';
        StageStitched_TextStringsArray(RowIndex, ColIndex).HandleToText = [];
        StageStitched_TextStringsArray(RowIndex, ColIndex).Color = [1 1 0];
    end
end



%START: PERFORM AUTOFOCUS AT OFFSET POSITION
AFCenterRowOffsetInMicrons = -GuiGlobalsStruct.MontageTarget.AF_Y_Offset_Microns;
AFCenterColOffsetInMicrons = GuiGlobalsStruct.MontageTarget.AF_X_Offset_Microns;

AFRowOffsetInMicrons = AFCenterRowOffsetInMicrons*r_target_north_UnitVector + ...
    AFCenterColOffsetInMicrons*c_target_north_UnitVector;
AFColOffsetInMicrons = AFCenterRowOffsetInMicrons*r_target_east_UnitVector +...
    AFCenterColOffsetInMicrons*c_target_east_UnitVector;

StageX_Meters = StageX_Meters_CenterOfMontage - AFColOffsetInMicrons/1000000;
StageY_Meters = StageY_Meters_CenterOfMontage - AFRowOffsetInMicrons/1000000;

MyStr = sprintf('Moving stage to(%0.5g, %0.5g)',StageX_Meters,StageY_Meters);
disp(MyStr);
GuiGlobalsStruct.MyCZEMAPIClass.MoveStage(StageX_Meters,StageY_Meters,stage_z,stage_t,stage_r,stage_m);
while(strcmp(GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeString('DP_STAGE_IS'),'Busy'))
    pause(.02)
end



if GuiGlobalsStruct.MontageParameters.IsSingle_AF_ForWholeMontage
    %reset original WD
    GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeSingle('AP_WD',GuiGlobalsStruct.WD_AtAcquitionStart);
    pause(1);
    %PerformAutoFocus;
    StartingMagForAF = GuiGlobalsStruct.MontageParameters.AutoFocusStartMag;
    IsPerformAutoStig = false;
    StartingMagForAS = round(StartingMagForAF/2);
    focOptions.IsDoQualCheck = GuiGlobalsStruct.MontageParameters.IsPerformQualityCheckOnEveryAF;
    focOptions.QualityThreshold = GuiGlobalsStruct.MontageParameters.AFQualityThreshold;
    Perform_AF_or_AFASAF(StartingMagForAF, IsPerformAutoStig, StartingMagForAS, focOptions);
end
if GuiGlobalsStruct.MontageParameters.IsSingle_AFASAF_ForWholeMontage
    %reset original WD + Stig
    GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeSingle('AP_WD',GuiGlobalsStruct.WD_AtAcquitionStart);
    BestGuess_StigX = median(GuiGlobalsStruct.StigX_ArrayOfValuesRecordedSinceStartOfMontageStack((max(1,end-(GuiGlobalsStruct.NumOfStigValuesToMedianOver-1)):end))); %takes median value of last 5 stigs
    BestGuess_StigY = median(GuiGlobalsStruct.StigY_ArrayOfValuesRecordedSinceStartOfMontageStack((max(1,end-(GuiGlobalsStruct.NumOfStigValuesToMedianOver-1)):end)));
    GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeSingle('AP_STIG_X',BestGuess_StigX);  %GuiGlobalsStruct.StigX_AtAcquitionStart);
    GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeSingle('AP_STIG_Y',BestGuess_StigY);  %GuiGlobalsStruct.StigY_AtAcquitionStart);
    pause(1);
    %PerformAutoFocusStigFocus;
    StartingMagForAF = GuiGlobalsStruct.MontageParameters.AutoFocusStartMag;
    IsPerformAutoStig = true;
    StartingMagForAS = round(StartingMagForAF/2);
    focOptions.IsDoQualCheck = GuiGlobalsStruct.MontageParameters.IsPerformQualityCheckOnEveryAF;
    focOptions.QualityThreshold = GuiGlobalsStruct.MontageParameters.AFQualityThreshold;
    Perform_AF_or_AFASAF(StartingMagForAF, IsPerformAutoStig, StartingMagForAS, focOptions);
    GuiGlobalsStruct.StigX_ArrayOfValuesRecordedSinceStartOfMontageStack(1+length(GuiGlobalsStruct.StigX_ArrayOfValuesRecordedSinceStartOfMontageStack)) = ...
        GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeSingle('AP_STIG_X'); %record this new value
    GuiGlobalsStruct.StigY_ArrayOfValuesRecordedSinceStartOfMontageStack(1+length(GuiGlobalsStruct.StigY_ArrayOfValuesRecordedSinceStartOfMontageStack)) = ...
        GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeSingle('AP_STIG_Y');
    
end
if GuiGlobalsStruct.MontageParameters.IsPlaneFit
    %reset original WD + Stig
    GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeSingle('AP_WD',GuiGlobalsStruct.WD_AtAcquitionStart);
    BestGuess_StigX = median(GuiGlobalsStruct.StigX_ArrayOfValuesRecordedSinceStartOfMontageStack((max(1,end-(GuiGlobalsStruct.NumOfStigValuesToMedianOver-1)):end))); %takes median value of last 5 stigs
    BestGuess_StigY = median(GuiGlobalsStruct.StigY_ArrayOfValuesRecordedSinceStartOfMontageStack((max(1,end-(GuiGlobalsStruct.NumOfStigValuesToMedianOver-1)):end)));
    GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeSingle('AP_STIG_X',BestGuess_StigX);  %GuiGlobalsStruct.StigX_AtAcquitionStart);
    GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeSingle('AP_STIG_Y',BestGuess_StigY);  %GuiGlobalsStruct.StigY_AtAcquitionStart);
    pause(1);
    RowDistanceBetweenTileCentersInMicrons_ForGridAutoFocus = GuiGlobalsStruct.MontageParameters.RowDistBetweenAFPointsMicrons; %50; %150;
    ColDistanceBetweenTileCentersInMicrons_ForGridAutoFocus = GuiGlobalsStruct.MontageParameters.ColDistBetweenAFPointsMicrons; %50; %150;
    ReturnedPlaneFitObject = GridAutoFocus_WithPlaneFit(RowDistanceBetweenTileCentersInMicrons_ForGridAutoFocus, ColDistanceBetweenTileCentersInMicrons_ForGridAutoFocus, MontageDirName);
    GuiGlobalsStruct.StigX_ArrayOfValuesRecordedSinceStartOfMontageStack(1+length(GuiGlobalsStruct.StigX_ArrayOfValuesRecordedSinceStartOfMontageStack)) = ...
        GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeSingle('AP_STIG_X'); %record this new value
    GuiGlobalsStruct.StigY_ArrayOfValuesRecordedSinceStartOfMontageStack(1+length(GuiGlobalsStruct.StigY_ArrayOfValuesRecordedSinceStartOfMontageStack)) = ...
        GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeSingle('AP_STIG_Y');
end



%Store the WD directly after this GridAutoFocus command to use as starting
%point for all others
StartingPointWD = GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeSingle('AP_WD');
StartingPoint_StigX = GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeSingle('AP_STIG_X');
StartingPoint_StigY = GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeSingle('AP_STIG_Y');


%Take test picture (remove when done testing)
% ImageWidthInPixels = 1000; %4096;%16384;%4096; %1024
% ImageHeightInPixels = 1000; %4096;%16384;%4096; %1024
% DwellTimeInMicroseconds = 2; %.5;
% FOV_microns = 100; %4nm pixel %GuiGlobalsStruct.MontageTarget.MontageWidthInMicrons; %40.96;
% IsDoAutoRetakeIfNeeded = false;
% IsMagOverride = false;
% MagForOverride = -1;
% WaferNameStr = WaferName;
% LabelStr = LabelStr;
%ImageFileNameStr = sprintf('%s\\AutoFocusPositionTestImage.tif', MontageDirName);
%Fibics_AcquireImage(ImageWidthInPixels, ImageHeightInPixels, DwellTimeInMicroseconds, ImageFileNameStr,...
%    FOV_microns, IsDoAutoRetakeIfNeeded, IsMagOverride, MagForOverride,  WaferNameStr, LabelStr);

%Added by AV on 01/07/2013. To blank the bean before moving to the center
%of Montage. This is to avaoud any beam damage.

GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeString('DP_BEAM_BLANKED','Yes'); % Added by AV 

StageX_Meters = StageX_Meters_CenterOfMontage;
StageY_Meters = StageY_Meters_CenterOfMontage;
MyStr = sprintf('Moving stage to(%0.5g, %0.5g)',StageX_Meters,StageY_Meters);
disp(MyStr);
GuiGlobalsStruct.MyCZEMAPIClass.MoveStage(StageX_Meters,StageY_Meters,stage_z,stage_t,stage_r,stage_m);
while(strcmp(GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeString('DP_STAGE_IS'),'Busy'))
    pause(.02)
end
%END: PERFORM AUTOFOCUS AT OFFSET POSITION

%pause (2);
%GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeString('DP_BEAM_BLANKED','No') % Added by AV
%pause (0.5);


%Take Montage Overview Image
IsTakeMontageOverviewImage = GuiGlobalsStruct.MontageParameters.IsAcquireOverviewImage;
if IsTakeMontageOverviewImage
    FOV_microns = GuiGlobalsStruct.MontageParameters.MontageOverviewImageFOV_microns;
    ImageWidthInPixels = GuiGlobalsStruct.MontageParameters.MontageOverviewImageWidth_pixels;
    ImageHeightInPixels = GuiGlobalsStruct.MontageParameters.MontageOverviewImageHeight_pixels;
    DwellTimeInMicroseconds = GuiGlobalsStruct.MontageParameters.MontageOverviewImageDwellTime_microseconds
    IsDoAutoRetakeIfNeeded = false;
    IsMagOverride = false;
    MagForOverride = -1;
    WaferNameStr = WaferName;
    LabelStr = LabelStr;
    
    ImageFileNameStr = sprintf('%s\\MontageOverviewImage_%s_sec%s.tif', MontageDirName, WaferName, LabelStr);
    Fibics_AcquireImage(ImageWidthInPixels, ImageHeightInPixels, DwellTimeInMicroseconds, ImageFileNameStr,...
        FOV_microns, IsDoAutoRetakeIfNeeded, IsMagOverride, MagForOverride,  WaferNameStr, LabelStr);
    
end


DropoutListFileName = sprintf('%s\\MontageTileDropOutList.txt',GuiGlobalsStruct.WaferDirectory);
if exist(DropoutListFileName,'file')
    DropOutListArray = dlmread(DropoutListFileName,',');
else
    DropOutListArray = [];
end

for RowIndex = NumRowTiles:1
    for ColIndex = NumColTiles:1
        
        IsDropOut = false;
        [NumDropOuts, dummy] = size(DropOutListArray);
        for DropOutListIndex = 1:NumDropOuts
            if (DropOutListArray(DropOutListIndex, 1) == RowIndex) && (DropOutListArray(DropOutListIndex, 2) == ColIndex)
                IsDropOut = true;
            end
        end
        
        if ~IsDropOut
            
            if (NumRowTiles == 1) && (NumColTiles == 1)
                ImageFileNameStr = sprintf('%s\\Image_%s.tif', MontageDirName, LabelStr);
            else
                %Tile_r1-c1_VMat2_sec01.tif
                ImageFileNameStr = sprintf('%s\\Tile_r%d-c%d_%s_sec%s.tif', MontageDirName, RowIndex, ColIndex, WaferName, LabelStr);
            end
            
            if ~exist(ImageFileNameStr, 'file') %do not take or move if already exists
                
                
                %%%% TAKE TILE IMAGE
                TileCenterRowOffsetInMicrons = (RowIndex -((NumRowTiles+1)/2)) * RowDistanceBetweenTileCentersInMicrons;
                TileCenterColOffsetInMicrons = (ColIndex -((NumColTiles+1)/2)) * ColDistanceBetweenTileCentersInMicrons;
                
                
                %Handle additional offset of full montage
                RowOffsetFromAlignTargetMicrons = -GuiGlobalsStruct.MontageTarget.YOffsetFromAlignTargetMicrons;
                ColOffsetFromAlignTargetMicrons = GuiGlobalsStruct.MontageTarget.XOffsetFromAlignTargetMicrons;
                TileCenterRowOffsetInMicrons = TileCenterRowOffsetInMicrons + RowOffsetFromAlignTargetMicrons;
                TileCenterColOffsetInMicrons = TileCenterColOffsetInMicrons + ColOffsetFromAlignTargetMicrons;
                %
                
                RowOffsetInMicrons = TileCenterRowOffsetInMicrons*r_target_north_UnitVector + ...
                    TileCenterColOffsetInMicrons*c_target_north_UnitVector;
                ColOffsetInMicrons = TileCenterRowOffsetInMicrons*r_target_east_UnitVector +...
                    TileCenterColOffsetInMicrons*c_target_east_UnitVector;
                
                
                
                StageX_Meters = StageX_Meters_CenterOfMontage - ColOffsetInMicrons/1000000;
                StageY_Meters = StageY_Meters_CenterOfMontage - RowOffsetInMicrons/1000000;
                
                MyStr = sprintf('Moving stage to(%0.5g, %0.5g)',StageX_Meters,StageY_Meters);
                disp(MyStr);
                GuiGlobalsStruct.MyCZEMAPIClass.MoveStage(StageX_Meters,StageY_Meters,stage_z,stage_t,stage_r,stage_m);
                while(strcmp(GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeString('DP_STAGE_IS'),'Busy'))
                    pause(.02)
                end
                
                pause(2);
                GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeString('DP_BEAM_BLANKED','No'); % Added by AV
                
                %Set the WD for this tile based on one of these three
                %options:
                if GuiGlobalsStruct.MontageParameters.IsAFOnEveryTile
                    GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeSingle('AP_WD',StartingPointWD);
                    pause(1); %1
                    %PerformAutoFocus;
                    StartingMagForAF = GuiGlobalsStruct.MontageParameters.AutoFocusStartMag;
                    IsPerformAutoStig = false;
                    StartingMagForAS = round(StartingMagForAF/2);
                    focOptions.IsDoQualCheck = GuiGlobalsStruct.MontageParameters.IsPerformQualityCheckOnEveryAF;
                    focOptions.QualityThreshold = GuiGlobalsStruct.MontageParameters.AFQualityThreshold;
                    Perform_AF_or_AFASAF(StartingMagForAF, IsPerformAutoStig, StartingMagForAS, focOptions);
                end
                

                if GuiGlobalsStruct.MontageParameters.IsAFASAFOnEveryTile
                    
                    if (RowIndex == NumRowTiles) && (ColIndex == NumColTiles)  % Added by AV 3/14/2013 to perform AF+AS+AS on the first tile only
                        
                        IsPerformAutoStig = true;
                        GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeSingle('AP_WD',StartingPointWD);
                        BestGuess_StigX = median(GuiGlobalsStruct.StigX_ArrayOfValuesRecordedSinceStartOfMontageStack((max(1,end-(GuiGlobalsStruct.NumOfStigValuesToMedianOver-1)):end))); %takes median value of last 5 stigs
                        BestGuess_StigY = median(GuiGlobalsStruct.StigY_ArrayOfValuesRecordedSinceStartOfMontageStack((max(1,end-(GuiGlobalsStruct.NumOfStigValuesToMedianOver-1)):end)));
                        GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeSingle('AP_STIG_X',BestGuess_StigX);  %GuiGlobalsStruct.StigX_AtAcquitionStart);
                        GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeSingle('AP_STIG_Y',BestGuess_StigY);  %GuiGlobalsStruct.StigY_AtAcquitionStart);
                        pause(1); %1
                        %PerformAutoFocusStigFocus;
                        StartingMagForAF = GuiGlobalsStruct.MontageParameters.AutoFocusStartMag;
                        StartingMagForAS = round(StartingMagForAF/2);
                        focOptions.IsDoQualCheck = GuiGlobalsStruct.MontageParameters.IsPerformQualityCheckOnEveryAF;
                        focOptions.QualityThreshold = GuiGlobalsStruct.MontageParameters.AFQualityThreshold;
                        Perform_AF_or_AFASAF(StartingMagForAF, IsPerformAutoStig, StartingMagForAS, focOptions);
                        GuiGlobalsStruct.StigX_ArrayOfValuesRecordedSinceStartOfMontageStack(1+length(GuiGlobalsStruct.StigX_ArrayOfValuesRecordedSinceStartOfMontageStack)) = ...
                            GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeSingle('AP_STIG_X'); %record this new value
                        GuiGlobalsStruct.StigY_ArrayOfValuesRecordedSinceStartOfMontageStack(1+length(GuiGlobalsStruct.StigY_ArrayOfValuesRecordedSinceStartOfMontageStack)) = ...
                            GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeSingle('AP_STIG_Y');
                    else
                        IsPerformAutoStig = false;
                        GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeSingle('AP_WD',StartingPointWD);
                        BestGuess_StigX = median(GuiGlobalsStruct.StigX_ArrayOfValuesRecordedSinceStartOfMontageStack((max(1,end-(GuiGlobalsStruct.NumOfStigValuesToMedianOver-1)):end))); %takes median value of last 5 stigs
                        BestGuess_StigY = median(GuiGlobalsStruct.StigY_ArrayOfValuesRecordedSinceStartOfMontageStack((max(1,end-(GuiGlobalsStruct.NumOfStigValuesToMedianOver-1)):end)));
                        %GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeSingle('AP_STIG_X',BestGuess_StigX);  %GuiGlobalsStruct.StigX_AtAcquitionStart);
                        %GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeSingle('AP_STIG_Y',BestGuess_StigY);  %GuiGlobalsStruct.StigY_AtAcquitionStart);
                        pause(1); %1
                        %PerformAutoFocusStigFocus;
                        StartingMagForAF = GuiGlobalsStruct.MontageParameters.AutoFocusStartMag;
                        StartingMagForAS = round(StartingMagForAF/2);
                        focOptions.IsDoQualCheck = GuiGlobalsStruct.MontageParameters.IsPerformQualityCheckOnEveryAF;
                        focOptions.QualityThreshold = GuiGlobalsStruct.MontageParameters.AFQualityThreshold;
                        Perform_AF_or_AFASAF(StartingMagForAF, IsPerformAutoStig, StartingMagForAS, focOptions);
                        GuiGlobalsStruct.StigX_ArrayOfValuesRecordedSinceStartOfMontageStack(1+length(GuiGlobalsStruct.StigX_ArrayOfValuesRecordedSinceStartOfMontageStack)) = ...
                            GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeSingle('AP_STIG_X'); %record this new value
                        GuiGlobalsStruct.StigY_ArrayOfValuesRecordedSinceStartOfMontageStack(1+length(GuiGlobalsStruct.StigY_ArrayOfValuesRecordedSinceStartOfMontageStack)) = ...
                            GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeSingle('AP_STIG_Y');
                    end
                end
     
                if GuiGlobalsStruct.MontageParameters.IsPlaneFit
                    if ~isempty(ReturnedPlaneFitObject)
                        NewWD = ReturnedPlaneFitObject(StageX_Meters,StageY_Meters);
                        %NewWD = ReturnedPlaneFitObject(StageX_Meters_CenterOfMontage,StageY_Meters_CenterOfMontage); %KH THIS ONLY DOWS AVERAGE AT CENTER!!!!!!!!!
                        GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeSingle('AP_WD',NewWD);
                        pause(1); %1
                    else
                        %%ErrorFileNameStr = sprintf('%s\\Error_PlaneFitReturnedEmptyMatrix_SkippingSection.mat', MontageDirName);
                        return; %KH added 11-14-2011 to skip entire section if plane fit failed
                    end
                end
                
                
                
                ImageWidthInPixels = GuiGlobalsStruct.MontageParameters.TileWidth_pixels; %4096;%16384;%4096; %1024
                ImageHeightInPixels = GuiGlobalsStruct.MontageParameters.TileWidth_pixels; %4096;%16384;%4096; %1024
                DwellTimeInMicroseconds = GuiGlobalsStruct.MontageParameters.TileDwellTime_microseconds; %.5;
                FOV_microns = GuiGlobalsStruct.MontageParameters.TileFOV_microns; %4nm pixel %GuiGlobalsStruct.MontageTarget.MontageWidthInMicrons; %40.96;
                IsDoAutoRetakeIfNeeded = false;
                %Fibics_AcquireImage(MyCZEMAPIClass, ImageWidthInPixels, ImageHeightInPixels, DwellTimeInMicroseconds, FileNameStr,...
                %      FOV_microns, IsDoAutoRetakeIfNeeded, IsMagOverride, MagForOverride,  WaferNameStr, LabelStr)
                IsMagOverride = false;
                MagForOverride = -1;
                WaferNameStr = WaferName;
                LabelStr = LabelStr;
                Fibics_AcquireImage(ImageWidthInPixels, ImageHeightInPixels, DwellTimeInMicroseconds, ImageFileNameStr,...
                    FOV_microns, IsDoAutoRetakeIfNeeded, IsMagOverride, MagForOverride,  WaferNameStr, LabelStr);
                
                %%Check Quality
                if exist('tilesTaken')
                    if exist('monQual')
                        startQ = length(monQual)+1;
                    else
                        startQ = 1;
                    end
                    for mQ = startQ : length(tilesTaken)
                        checkTile = tilesTaken{mQ};
                        qual = checkFileQual(checkTile);
                        monQual(mQ) = qual.quality
                        q.fileNum = q.fileNum+1;
                        q.qual(q.fileNum) = qual;
                        q.tile{q.fileNum} = checkTile;
                        q.s(str2num(LabelStr)).r(RowIndex).c(ColIndex) = qual; %KH Ask JM if this R,C is right???
                        save([GuiGlobalsStruct.TempImagesDirectory '\watchQ.mat'],'q')
                        
                        %KH Update StageStitched text labels                    
                        StageStitched_TextStringsArray(tilesTaken_RowNum(mQ), tilesTaken_ColNum(mQ)).Text = ...
                            sprintf('%0.4g',qual.quality);
                        if qual.quality >= GuiGlobalsStruct.MontageParameters.ImageQualityThreshold
                            StageStitched_TextStringsArray(tilesTaken_RowNum(mQ), tilesTaken_ColNum(mQ)).Color = [0 1 0];
                        else
                            StageStitched_TextStringsArray(tilesTaken_RowNum(mQ), tilesTaken_ColNum(mQ)).Color = [1 0 0];
                        end
                    end
                end
                
                %%Wait for image to be finished
                while(GuiGlobalsStruct.MyCZEMAPIClass.Fibics_IsBusy)
                    pause(.2); %1
                end
                
                
                %Try to read this image downsampled. Continue trying until it works meaning
                %that Fibics is done writing file
                IfCheckFileIsWritten = true;
                if IfCheckFileIsWritten
                    
                    IsReadOK = false;
                    while ~IsReadOK
                        IsReadOK = true;
                        try
                            MyDownSampledImage = imread(ImageFileNameStr, 'PixelRegion', {[1 16 ImageWidthInPixels], [1 16 ImageWidthInPixels]});
                        catch MyException
                            IsReadOK = false;
                            pause(0.5);
                        end
                    end
                    
                    %%Done waiting
                    MyNewIndex = length(tilesTaken) + 1;
                    tilesTaken{MyNewIndex} = ImageFileNameStr;% Add new image to list
                    tilesTaken_RowNum(MyNewIndex) = RowIndex;
                    tilesTaken_ColNum(MyNewIndex) = ColIndex;
                    
                    IsDisplay = true;
                    if IsDisplay
                        %imread(ImageFileNameStr, 'tif', 'PixelRegion',{[START INCREMENT STOP], [START INCREMENT STOP]});
                        MyImage = imread(ImageFileNameStr, 'tif', 'PixelRegion', {[1 Increment ImageHeightInPixels],[1 Increment ImageWidthInPixels]});
                        
                        %put in border. Remember this image is just for show, it
                        %does not even compensate for the tile overlaps
                        MyImage(1:BorderPixels,:) = 0;
                        MyImage(end-BorderPixels+1:end,:) = 0;
                        MyImage(:,1:BorderPixels) = 0;
                        MyImage(:,end-BorderPixels+1:end) = 0;
                        
                        [MaxSSTileR, MaxSSTileC] = size(MyImage);
                        StartR = (MaxSSTileR*(RowIndex-1))+1;
                        StartC = (MaxSSTileC*(ColIndex-1))+1;
                        StageStitchedImage(StartR:StartR+MaxSSTileR-1, StartC:StartC+MaxSSTileC-1) = MyImage;
                        figure(StitchFigNum);
                        clf;
                        title(MontageDirName);
                        %subplot(NumRowTiles, NumColTiles, ColIndex + NumColTiles*(RowIndex-1));
                        imshow(StageStitchedImage,[0, 255]);
                        
                        
                        StageStitched_TextStringsArray(RowIndex, ColIndex).textX = StartC+(MaxSSTileC/2);
                        StageStitched_TextStringsArray(RowIndex, ColIndex).textY = StartR+(MaxSSTileR/2);
                        StageStitched_TextStringsArray(RowIndex, ColIndex).Text = sprintf('(%d, %d)',RowIndex, ColIndex);
                        
                        UpdateTextOnStageStitched(NumRowTiles, NumColTiles, StitchFigNum, StageStitched_TextStringsArray);
                        
                    end
                end
            end
            %         r_target_offset = r_target + TileCenterRowOffsetInPixels*r_target_north_UnitVector + TileCenterColOffsetInPixels*c_target_north_UnitVector;
            %         c_target_offset = c_target + TileCenterRowOffsetInPixels*r_target_east_UnitVector + TileCenterColOffsetInPixels*c_target_east_UnitVector;
            
        end
    end
end

%check quality of last  if exist('tilesTaken')
if exist('tilesTaken')
    if exist('monQual')
        startQ = length(monQual)+1;
    else
        startQ = 1;
    end
    for mQ = startQ : length(tilesTaken)
        checkTile = tilesTaken{mQ};
        qual = checkFileQual(checkTile);
        monQual(mQ) = qual.quality
        q.fileNum = q.fileNum+1;
        qual
        q
        q.qual(q.fileNum) = qual;
        q.tile{q.fileNum} = checkTile;
        save([GuiGlobalsStruct.TempImagesDirectory '\watchQ.mat'],'q')
        
        %KH Update StageStitched text labels
        StageStitched_TextStringsArray(tilesTaken_RowNum(mQ), tilesTaken_ColNum(mQ)).Text = ...
            sprintf('%0.4g',qual.quality);
        if qual.quality >= GuiGlobalsStruct.MontageParameters.ImageQualityThreshold
            StageStitched_TextStringsArray(tilesTaken_RowNum(mQ), tilesTaken_ColNum(mQ)).Color = [0 1 0];
        else
            StageStitched_TextStringsArray(tilesTaken_RowNum(mQ), tilesTaken_ColNum(mQ)).Color = [1 0 0];
        end
    end
end

%do update of ss display
figure(StitchFigNum);
clf;
imshow(StageStitchedImage,[0, 255]);
UpdateTextOnStageStitched(NumRowTiles, NumColTiles, StitchFigNum, StageStitched_TextStringsArray);




%%RETAKE BAD IMAGES
if GuiGlobalsStruct.MontageParameters.IsPerformQualCheckAfterEachImage == true
    
    %%Find bad files
    retakeTiles = {tilesTaken{find(monQual <= GuiGlobalsStruct.MontageParameters.ImageQualityThreshold)}}
    
    if length(retakeTiles) > 0
        disp('About to retake the following tiles that did not pass the qual check:');
        for RetakeNum = 1:length(retakeTiles)
            disp(sprintf('   %s', retakeTiles{RetakeNum}));
        end
        disp(' ');
        
        
        for RetakeNum = 1:length(retakeTiles)
            ImageFileNameStr = retakeTiles{RetakeNum};
            
            
            
            
            
            
            DataFileNameStr = sprintf('%s.mat', ImageFileNameStr(1:length(ImageFileNameStr)-4));
            disp(sprintf('Retaking:   %s', ImageFileNameStr));
            disp(sprintf('   Data file name:   %s', DataFileNameStr));
            
            %load in data file with position
            load(DataFileNameStr, 'Info');
            StageX_Meters = Info.StageX_Meters;
            StageY_Meters = Info.StageY_Meters;
            
            %move to tile position
            MyStr = sprintf('Moving stage to(%0.5g, %0.5g)',StageX_Meters,StageY_Meters);
            disp(MyStr);
            GuiGlobalsStruct.MyCZEMAPIClass.MoveStage(StageX_Meters,StageY_Meters,stage_z,stage_t,stage_r,stage_m);
            while(strcmp(GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeString('DP_STAGE_IS'),'Busy'))
                pause(.02)
            end
            
            %perform autofocus
            GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeSingle('AP_STIG_X',BestGuess_StigX);
            GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeSingle('AP_STIG_Y',BestGuess_StigY); 
            GuiGlobalsStruct.MyCZEMAPIClass.Set_PassedTypeSingle('AP_WD',StartingPointWD);
            pause(1); %1
            StartingMagForAF = GuiGlobalsStruct.MontageParameters.AutoFocusStartMag;
            IsPerformAutoStig = true; %false;
            StartingMagForAS = StartingMagForAF; %round(StartingMagForAF/2);
            focOptions.IsDoQualCheck = GuiGlobalsStruct.MontageParameters.IsPerformQualityCheckOnEveryAF;
            focOptions.QualityThreshold = GuiGlobalsStruct.MontageParameters.AFQualityThreshold;
            Perform_AF_or_AFASAF(StartingMagForAF, IsPerformAutoStig, StartingMagForAS, focOptions);
            
            %take image
            ImageWidthInPixels = GuiGlobalsStruct.MontageParameters.TileWidth_pixels; %4096;%16384;%4096; %1024
            ImageHeightInPixels = GuiGlobalsStruct.MontageParameters.TileWidth_pixels; %4096;%16384;%4096; %1024
            DwellTimeInMicroseconds = GuiGlobalsStruct.MontageParameters.TileDwellTime_microseconds; %.5;
            FOV_microns = GuiGlobalsStruct.MontageParameters.TileFOV_microns; %4nm pixel %GuiGlobalsStruct.MontageTarget.MontageWidthInMicrons; %40.96;
            IsDoAutoRetakeIfNeeded = false;
            IsMagOverride = false;
            MagForOverride = -1;
            WaferNameStr = WaferName;
            LabelStr = LabelStr;
            Fibics_AcquireImage(ImageWidthInPixels, ImageHeightInPixels, DwellTimeInMicroseconds, ImageFileNameStr,...
                FOV_microns, IsDoAutoRetakeIfNeeded, IsMagOverride, MagForOverride,  WaferNameStr, LabelStr);
            %%Wait for image to be finished
            while(GuiGlobalsStruct.MyCZEMAPIClass.Fibics_IsBusy)
                pause(.2); %1
            end
            %Try to read this image downsampled. Continue trying until it works meaning
            %that Fibics is done writing file
            IfCheckFileIsWritten = true;
            if IfCheckFileIsWritten
                IsReadOK = false;
                while ~IsReadOK
                    IsReadOK = true;
                    try
                        MyDownSampledImage = imread(ImageFileNameStr, 'PixelRegion', {[1 16 ImageWidthInPixels], [1 16 ImageWidthInPixels]});
                    catch MyException
                        IsReadOK = false;
                        pause(0.5);
                    end
                end
            end
            
            
            
            %Incorporate into stage stitched and display
            if IsDisplay
                %imread(ImageFileNameStr, 'tif', 'PixelRegion',{[START INCREMENT STOP], [START INCREMENT STOP]});
                MyImage = imread(ImageFileNameStr, 'tif', 'PixelRegion', {[1 Increment ImageHeightInPixels],[1 Increment ImageWidthInPixels]});
                
                %put in border. Remember this image is just for show, it
                %does not even compensate for the tile overlaps
                MyImage(1:BorderPixels,:) = 0;
                MyImage(end-BorderPixels+1:end,:) = 0;
                MyImage(:,1:BorderPixels) = 0;
                MyImage(:,end-BorderPixels+1:end) = 0;
                
                %Extract row anbd col numbers from file name
                A1 = findstr(ImageFileNameStr, 'Tile_r');
                RowStartIndex = A1(end)+ 6;
                A2 = findstr(ImageFileNameStr(RowStartIndex:end), '-c');
                RowEndIndex = RowStartIndex + A2(1) - 2;
                A3 = findstr(ImageFileNameStr, '-c');
                ColStartIndex = A3(end) + 2;
                A4 = findstr(ImageFileNameStr(ColStartIndex:end), '_');
                ColEndIndex = ColStartIndex + A4(1) - 2;
                RowStr = ImageFileNameStr(RowStartIndex:RowEndIndex);
                ColStr = ImageFileNameStr(ColStartIndex:ColEndIndex);
                RowIndex = str2num(RowStr);
                ColIndex = str2num(ColStr);
                
                
                [MaxSSTileR, MaxSSTileC] = size(MyImage);
                StartR = (MaxSSTileR*(RowIndex-1))+1;
                StartC = (MaxSSTileC*(ColIndex-1))+1;
                StageStitchedImage(StartR:StartR+MaxSSTileR-1, StartC:StartC+MaxSSTileC-1) = MyImage;
                figure(StitchFigNum);
                clf;
                title(MontageDirName);
                %subplot(NumRowTiles, NumColTiles, ColIndex + NumColTiles*(RowIndex-1));
                imshow(StageStitchedImage,[0, 255]);
                
                
                StageStitched_TextStringsArray(RowIndex, ColIndex).textX = StartC+(MaxSSTileC/2);
                StageStitched_TextStringsArray(RowIndex, ColIndex).textY = StartR+(MaxSSTileR/2);
                StageStitched_TextStringsArray(RowIndex, ColIndex).Text = sprintf('(%d, %d)',RowIndex, ColIndex);
                
                %check quality immediatly
                qual = checkFileQual(ImageFileNameStr);
                
                %KH Update StageStitched text labels
                StageStitched_TextStringsArray(RowIndex, ColIndex).Text = ...
                    sprintf('%0.4g',qual.quality);
                if qual.quality >= GuiGlobalsStruct.MontageParameters.ImageQualityThreshold
                    StageStitched_TextStringsArray(RowIndex, ColIndex).Color = [0 1 0];
                else
                    StageStitched_TextStringsArray(RowIndex, ColIndex).Color = [1 0 0];
                end
                
                UpdateTextOnStageStitched(NumRowTiles, NumColTiles, StitchFigNum, StageStitched_TextStringsArray);
                
            end
            
        end
        
        
    end
    disp(' ');
    
    
else
    disp('All tiles on this section passed the qual check. No retakes needed.');
end


StageStitchedImageFileNameStr = sprintf('%s\\StageStitched_%s_sec%s.tif', MontageDirName, WaferName, LabelStr);
imwrite(StageStitchedImage, StageStitchedImageFileNameStr, 'tif');

if GuiGlobalsStruct.MontageParameters.IsPerformQualCheckAfterEachImage == true
    StageStitchedImageWithQualValsFileNameStr = sprintf('%s\\StageStitched_%s_sec%s_WithQualVals.tif', MontageDirName, WaferName, LabelStr);
    saveas(StitchFigNum,StageStitchedImageWithQualValsFileNameStr,'tif');
end

%MOVE BACK TO ORIGINAL POSITION
StageX_Meters = StageX_Meters_CenterOfMontage;
StageY_Meters = StageY_Meters_CenterOfMontage;

MyStr = sprintf('Moving stage to(%0.5g, %0.5g)',StageX_Meters,StageY_Meters);
disp(MyStr);
GuiGlobalsStruct.MyCZEMAPIClass.MoveStage(StageX_Meters,StageY_Meters,stage_z,stage_t,stage_r,stage_m);
while(strcmp(GuiGlobalsStruct.MyCZEMAPIClass.Get_ReturnTypeString('DP_STAGE_IS'),'Busy'))
    pause(.02)
end



end

function UpdateTextOnStageStitched(NumRowTiles, NumColTiles, StitchFigNum, StageStitched_TextStringsArray)

for tR = 1:NumRowTiles
    for tC = 1:NumColTiles
        figure(StitchFigNum);
        
        textX = StageStitched_TextStringsArray(tR, tC).textX;
        textY = StageStitched_TextStringsArray(tR, tC).textY;
        MyText = StageStitched_TextStringsArray(tR, tC).Text;
        %TEXT(X,Y,'string')
        h = text(textX, textY, MyText);
        set(h,'Color', StageStitched_TextStringsArray(tR, tC).Color);
        set(h,'FontSize', 18);
    end
end

end

    