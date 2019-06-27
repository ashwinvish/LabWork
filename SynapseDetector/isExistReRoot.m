function [logical] = isExistReRoot(cellID)
% isExistReRoot checks if cellID has be reRooted 

filePath  = '/Users/ashwin/Google Drive/Zfish/LowEMtoHighEM/SWC_all/consensus-20180920/swc/';
reSampleFactor = 5000;

for i = 1:length(cellID)
    fileName = sprintf('%d_reRoot_reSample_%d.swc',cellID(i),reSampleFactor);
    if exist (fullfile(filePath, fileName))
        logical(i) = true;
    else
        logical(i) = false;
    end
end

end