clc;
clear all;

MontageDir = uigetdir('.','Pick the directiry');
if MontageDir == 0
    return;
end

pwd = MontageDir;

for r = 1:4
    for c = 1:4
        Fname = 'Tile_r%d-c%d_S2-W008_sec1.mat';
        str = sprintf(Fname,r,c);
        Param(r,c) = load(fullfile(MontageDir, str));
        X(r,c) = Param(r,c).Info.StageX_Meters;
        Y(r,c) = Param(r,c).Info.StageY_Meters;
        Z(r,c) = Param(r,c).Info.WorkingDistance;
        
    end
end

MedianValueGridAF_WD = median(Z);

ModelStr = 'poly11'; %'poly11'
%generate plane fit object
ReturnedPlaneFitObject = fit( [X, Y], ...
    Z, ModelStr); %'lowess'