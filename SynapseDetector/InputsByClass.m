function [cellClasses] = InputsByClass(cellID,df);

[cellClasses.Inputs, cellClasses.PSDID] = SynapticPartners(cellID,1,df);

for i =1:size(cellClasses.PSDID,1)
    cellClasses.PSDsize(i,1) = df.size(df.psd_segid == cellClasses.PSDID(i));
end

cellClasses.PreSynCoords = PrePartnerCoordinates(cellClasses.PSDID,df);

cellClasses.isSaccadic = isSaccade(cellClasses.Inputs);
cellClasses.isVestibular = isVestibular(cellClasses.Inputs);
cellClasses.isContra = isContra(cellClasses.Inputs);
cellClasses.isIntegrator = isIntegrator(cellClasses.Inputs);

cellClasses.Saccadic = cellClasses.Inputs(cellClasses.isSaccadic);
cellClasses.Vestibular = cellClasses.Inputs(cellClasses.isVestibular);
cellClasses.Contra  = cellClasses.Inputs(cellClasses.isContra);
cellClasses.Integrator = cellClasses.Inputs(cellClasses.isIntegrator);

cellClasses.EverythingElse = setdiff(cellClasses.Inputs, ...
    [cellClasses.Saccadic;cellClasses.Vestibular;cellClasses.Contra;cellClasses.Integrator]);
cellClasses.isEverythingElse = ismember(cellClasses.Inputs,cellClasses.EverythingElse);

cellClasses.MotorDist = isMotor(cellID,df);

fname  = '/Users/ashwin/Documents/LowEMtoHighEM/SWC_all/consensus-20180920/swc/';
filename = sprintf('%d.swc',cellID);

if exist(fullfile(fname,filename))>0
    temp = dlmread(fullfile(fname,filename), ' ');
    coord = temp(:,3:5);
    cellClasses.Origin = TransformPoints(coord(1,:),4);
else
    cellClasses.Origin = NaN;
end



end
