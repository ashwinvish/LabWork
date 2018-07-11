fname = '/usr/people/ashwinv/seungmount/research/Ashwin/Z-Brain/LowEMtoHighEM/SWC_all/consensus-20180412/';
cd(fname);
directory = dir('*.swc');
for i = 1:size(directory,1)
    [ContraCells(i).isContra,ContraCells(i).AngleinDegrees] = isMidline(directory(i).name(1:5),false);
    ContraCells(i).CellID = directory(i).name(1:5);
end
contra = [];
ipsi = [];
for i = 1:size(ContraCells,2)
    if ContraCells(i).isContra==1
        contra = [contra; str2num(ContraCells(i).CellID)];
    else
        ipsi = [ipsi;str2num(ContraCells(i).CellID)];
    end
end