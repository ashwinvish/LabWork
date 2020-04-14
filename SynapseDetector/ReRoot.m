function[] = ReRoot(cellID)
% re-root trees

 %if isExistReRoot(cellID) == 1
 %    disp('Rerooted file exists');
 %else
    
    filePath  = '/Users/ashwin/Google Drive/Zfish/LowEMtoHighEM/SWC_all/consensus-20190415/swc/';
    oldFilePath = '/Users/ashwin/Google Drive/Zfish/LowEMtoHighEM/SWC_all/combinedConsensus-resampled/';
    fileName = sprintf('%d.swc',cellID);
    
    if exist(fullfile(filePath,fileName))==0
        filePath  = '/Users/ashwin/Google Drive/Zfish/LowEMtoHighEM/SWC_all/consensus-20190923/swc/';
        oldFilePath = '/Users/ashwin/Google Drive/Zfish/LowEMtoHighEM/SWC_all/combinedConsensus-resampled/';
    end
    
    reSampleFactor = 5000;
    
    if exist(fullfile(filePath,fileName))
        [tree,~,~] = load_tree(fullfile(filePath,fileName));
        
        tree = resample_tree(tree,reSampleFactor);
        xplore_tree(tree);
        view(90,0);
        prompt = 'Enter the new root node: '
        %lia = find(tree.X >= max(tree.X)-10000,1);    
        %pmt = find(tree.Y == min(tree.Y));
        %pmt = 1;
        
        
        %if ~isempty(pmt) 
        
        
        istart = input(prompt);
        %istart = pmt;
        close all;
        [tree, order] = redirect_tree (tree, istart);
        newFileName = sprintf('%d_reRoot_reSample_%d.swc',cellID,reSampleFactor);
        swc_tree (tree, fullfile(oldFilePath,newFileName));
    else
        close all;
        newFileName = sprintf('%d_ORoot_reSample_%d.swc',cellID,reSampleFactor);
        swc_tree (tree, fullfile(oldFilePath,newFileName));
    end
    
        
    %end
%end

end
