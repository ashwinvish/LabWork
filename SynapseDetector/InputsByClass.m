function [cellProperties] = InputsByClass(cellID,df);

% transfrom teee to Z-brian space, in microns
cellProperties.Tree = SwctoZbrian(cellID);
cellProperties.Origin = [cellProperties.Tree{1}.X(1), cellProperties.Tree{1}.Y(1), cellProperties.Tree{1}.Z(1)];
cellProperties.Rhombomere = struct2array(isRhombomere(cellID));
[cellProperties.Inputs, cellProperties.PSDID] = SynapticPartners(cellID,1,df);
%cellProperties.InputsRhombomeres = struct2array(isRhombomere(cellProperties.Inputs));

for i =1:size(cellProperties.PSDID,1)
    cellProperties.PSDsize(i,1) = df.size(df.psd_segid == cellProperties.PSDID(i));
end

cellProperties.PreSynCoords = PrePartnerCoordinates(cellProperties.PSDID,df);
cellProperties.PreSynCoordsTransformed = TransformPoints(cellProperties.PreSynCoords,0);
cellProperties.PathLength =  PathLengthToCoordinate(cellProperties.PreSynCoordsTransformed,cellProperties.Tree{1});

cellProperties.isSaccadic = isSaccade(cellProperties.Inputs);
cellProperties.isVestibular = isVestibular(cellProperties.Inputs);
cellProperties.isContra = isContra(cellProperties.Inputs);
cellProperties.isIntegrator = isIntegrator(cellProperties.Inputs);

cellProperties.Saccadic = cellProperties.Inputs(cellProperties.isSaccadic);
cellProperties.Vestibular = cellProperties.Inputs(cellProperties.isVestibular);
cellProperties.Contra  = cellProperties.Inputs(cellProperties.isContra);
cellProperties.Integrator = cellProperties.Inputs(cellProperties.isIntegrator);


idx = ~ismember(cellProperties.Inputs, ...
    [cellProperties.Saccadic;cellProperties.Vestibular;cellProperties.Contra;cellProperties.Integrator],'rows');
cellProperties.EverythingElse = cellProperties.Inputs(idx);
cellProperties.isEverythingElse = ismember(cellProperties.Inputs,cellProperties.EverythingElse);

cellProperties.MotorDist = isMotor(cellID,df);

end
