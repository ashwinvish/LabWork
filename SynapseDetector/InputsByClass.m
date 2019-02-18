function [cellClasses] = InputsByClass(cellID,df);

% transfrom teee to Z-brian space, in microns
cellClasses.Tree = SwctoZbrian(cellID);
cellClasses.Origin = [cellClasses.Tree{1}.X(1), cellClasses.Tree{1}.Y(1), cellClasses.Tree{1}.Z(1)];
[cellClasses.Inputs, cellClasses.PSDID] = SynapticPartners(cellID,1,df);

for i =1:size(cellClasses.PSDID,1)
    cellClasses.PSDsize(i,1) = df.size(df.psd_segid == cellClasses.PSDID(i));
end

cellClasses.PreSynCoords = PrePartnerCoordinates(cellClasses.PSDID,df);
cellClasses.PreSynCoordsTransformed = TransformPoints(cellClasses.PreSynCoords,0);
cellClasses.PathLength =  PathLengthToCoordinate(cellClasses.PreSynCoordsTransformed,cellClasses.Tree{1});

cellClasses.isSaccadic = isSaccade(cellClasses.Inputs);
cellClasses.isVestibular = isVestibular(cellClasses.Inputs);
cellClasses.isContra = isContra(cellClasses.Inputs);
cellClasses.isIntegrator = isIntegrator(cellClasses.Inputs);

cellClasses.Saccadic = cellClasses.Inputs(cellClasses.isSaccadic);
cellClasses.Vestibular = cellClasses.Inputs(cellClasses.isVestibular);
cellClasses.Contra  = cellClasses.Inputs(cellClasses.isContra);
cellClasses.Integrator = cellClasses.Inputs(cellClasses.isIntegrator);

idx = ~ismember(cellClasses.Inputs, ...
    [cellClasses.Saccadic;cellClasses.Vestibular;cellClasses.Contra;cellClasses.Integrator],'rows');
cellClasses.EverythingElse = cellClasses.Inputs(idx);
cellClasses.isEverythingElse = ismember(cellClasses.Inputs,cellClasses.EverythingElse);

cellClasses.MotorDist = isMotor(cellID,df);

end
