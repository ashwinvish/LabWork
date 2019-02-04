function[] = ReRoot(cellID)
% re-root trees

filePath  = '/Users/ashwin/Documents/LowEMtoHighEM/SWC_all/consensus-20180920/swc/';
fileName = sprintf('%d.swc',cellID);

[tree,~,~] = load_tree(fullfile(filePath,fileName));
reSampleFactor = 5000;
tree = resample_tree(tree,reSampleFactor);
xplore_tree(tree);
prompt = 'Enter the new root node: '

if ~isempty(prompt)
    istart = input(prompt);
    close all;
    [tree, order] = redirect_tree (tree, istart);
    newFileName = sprintf('%d_reRoot_reSample_%d.swc',cellID,reSampleFactor);
    swc_tree (tree, fullfile(filePath,newFileName));
else
    close all;
    newFileName = sprintf('%d_ORoot_reSample_%d.swc',cellID,reSampleFactor);
    swc_tree (tree, fullfile(filePath,newFileName));
end

end