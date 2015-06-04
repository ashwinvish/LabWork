   global GuiGlobalsStruct
    Directory = uigetdir('Z:\MasterUTSLdirectory\Ashwin\');
    FOV_microns = 12.05; % in microns
    ImageWidthInPixels = 4000;
    ImageHeightInPixels = 4000;
    DwellTimeInMicroseconds = [0.1;0.2;0.5;1;1.5;2.0];
    IsDoAutoRetakeIfNeeded = false;
    IsMagOverride = false;
    MagForOverride = -1;
    WaferNameStr = 'S2-W001';
    LabelStr = '141';
    EHT = '5kV';
    WD = '6600'; % in microns
    imageno = 1;
    %clear imageno ;
   
    
    
    
    MontageDirName = sprintf('%s\\%s', Directory,EHT);
    
    
     disp(sprintf('Creating directory: %s',MontageDirName));
     [success,message,messageid] = mkdir(MontageDirName);
     
     %imageno = length(DwellTimeInMicroseconds);
     
     for imageno = 1:length(DwellTimeInMicroseconds)
         
         disp(sprintf('Creating File: %s\\%03d_SNROverviewImage_%s_%s_%s.tif', MontageDirName,imageno, EHT, WD ,num2str(DwellTimeInMicroseconds(imageno)*1000 )));
         ImageFileNameStr = sprintf('%s\\%03d_SNROverviewImage_%s_%s_%s.tif', MontageDirName,imageno, EHT, WD ,num2str(DwellTimeInMicroseconds(imageno)*1000));
         Fibics_AcquireImage(ImageWidthInPixels, ImageHeightInPixels, DwellTimeInMicroseconds(imageno), ImageFileNameStr,...
             FOV_microns, IsDoAutoRetakeIfNeeded, IsMagOverride, MagForOverride, WaferNameStr, LabelStr);
         
         pause(1);
         %imageno = length(DwellTimeInMicroseconds) - temp ;
     end