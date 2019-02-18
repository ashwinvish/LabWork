function[] = ReRoot(cellID)
% re-root trees

filePath  = '/Users/ashwin/Documents/LowEMtoHighEM/SWC_all/consensus-20180920/swc/';
oldFilePath = filePath;
fileName = sprintf('%d.swc',cellID);

if exist(fullfile(filePath,fileName))==0
    filePath  = '/Users/ashwin/Documents/LowEMtoHighEM/SWC_all/consensus-20181125/';
    oldFilePath = '/Users/ashwin/Documents/LowEMtoHighEM/SWC_all/consensus-20180920/swc/';
end

if exist(fullfile(filePath,fileName)) 
    [tree,~,~] = load_tree(fullfile(filePath,fileName));
    reSampleFactor = 5000;
    tree = resample_tree(tree,reSampleFactor);
    xplore_tree(tree);
    view(90,0);
    prompt = 'Enter the new root node: '
    
    %lia = find(tree.X >= max(tree.X)-10000,1);
    %pmt = find(tree.Y == min(tree.Y));
    
    if ~isempty(prompt)
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

    
end

end
