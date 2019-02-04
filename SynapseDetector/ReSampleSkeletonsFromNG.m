clc;
clear;
if ismac
    addpath(genpath('/Users/ashwin/Documents/LabWork'));
    filePath  = '/Users/ashwin/Documents/LowEMtoHighEM/SWC_all/consensus-20180920/swc/';
else
    addpath(genpath('/usr/people/ashwinv/seungmount/research/Ashwin/Scripts'));
    filepath = '/usr/people/ashwinv/seungmount/research/Ashwin/Z-Brain/LowEMtoHighEM/SWC_all/consensus-20180920/swc/'
end
    cd(filePath);
    D = dir('*.swc');
    
    parfor i = 1:size(D,1)
        if (strfind(D(i).name,'.') == 6)
            fileName = sprintf('%s',D(i).name);
            fprintf('processing files : %s', fileName)
            [tree,~,~] = load_tree(fullfile(filePath,fileName));
            reSampleFactor = 5000;
            tree = resample_tree(tree,reSampleFactor);
            newFileName = sprintf('%s_reSample_%d.swc',D(i).name,reSampleFactor);
            swc_tree(tree, fullfile(filePath,newFileName));t
        else
            continue;
        end
        
    end