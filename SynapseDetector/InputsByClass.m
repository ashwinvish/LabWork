function [cellProperties] = InputsByClass(cellID,df,synapseIden);
 %synapseIden = 2; % look for postsynaptic neurons
 %synapseIden = 1; % presynaptic neurons, default

% if nargin >2
%     synapseIden = 2; % look for postsynaptic neurons
% else
%     synapseIden = 1; % presynaptic neurons, default
% end

cellProperties.cellID = cellID;


% transfrom teee to Z-brian space, in microns
if (isExistReRoot(cellID) == 1)
    cellProperties.Tree = SwctoZbrian(cellID);
    cellProperties.Origin = [cellProperties.Tree{1}.X(1), cellProperties.Tree{1}.Y(1), cellProperties.Tree{1}.Z(1)];
    cellProperties.Rhombomere = struct2array(isRhombomere(cellID));
else
    cellProperties.Tree = [];
    cellProperties.Origin = [];
    cellProperties.Rhombomere = [];
end


if synapseIden == 2
    [cellProperties.Outputs, cellProperties.PSDID] = SynapticPartners(cellID,synapseIden,df);
    if ~isempty(cellProperties.Outputs)
        for i =1:size(cellProperties.PSDID,1)
            cellProperties.PSDsize(i,1) = df.size(df.psd_segid == cellProperties.PSDID(i));
        end
        
        cellProperties.PostSynCoords = PrePartnerCoordinates(cellProperties.PSDID,df);
        cellProperties.PostSynCoordsTransformed = TransformPoints(cellProperties.PostSynCoords,0);
        
        if (isExistReRoot(cellID) ==1 )
            cellProperties.PathLength =  PathLengthToCoordinate(cellProperties.PostSynCoordsTransformed,cellProperties.Tree{1});
        else
            cellProperties.PathLength = [];
        end
        
        cellProperties.isSaccadic = isSaccade(cellProperties.Outputs);
        cellProperties.isVestibular = isVestibular(cellProperties.Outputs);
        cellProperties.isContra = isContra(cellProperties.Outputs);
        cellProperties.isIntegrator = isIntegrator(cellProperties.Outputs);
        
        cellProperties.Saccadic = cellProperties.Outputs(cellProperties.isSaccadic);
        cellProperties.Vestibular = cellProperties.Outputs(cellProperties.isVestibular);
        cellProperties.Contra  = cellProperties.Outputs(cellProperties.isContra);
        cellProperties.Integrator = cellProperties.Outputs(cellProperties.isIntegrator);
        
        
%         idx = ~ismember(cellProperties.Outputs, ...
%             [cellProperties.Saccadic;cellProperties.Vestibular;cellProperties.Contra;cellProperties.Integrator],'rows');
%        cellProperties.EverythingElse = cellProperties.Outputs(idx);

        cellProperties.EverythingElse  = setdiff(cellProperties.Outputs,[cellProperties.Saccadic;cellProperties.Vestibular;cellProperties.Contra;cellProperties.Integrator]);
        cellProperties.isEverythingElse = ismember(cellProperties.Outputs,cellProperties.EverythingElse);
        
        cellProperties.MotorDist = isMotor(cellID,df); % 1X5 (cellID, ABDr, ABDc, ABDir, ABDic)
    else
        cellProperties.PSDsize = [];
        cellProperties.PostSynCoords =[];
        cellProperties.PostSynCoordsTransformed=[];
        cellProperties.PathLength =[];
        cellProperties.isSaccadic =[];
        cellProperties.isVestibular = [];
        cellProperties.isContra = [];
        cellProperties.isIntegrator = [];
        cellProperties.Saccadic = [];
        cellProperties.Vestibular = [];
        cellProperties.Contra  = [];
        cellProperties.Integrator = [];
        cellProperties.EverythingElse =[];
        cellProperties.isEverythingElse = [];
        cellProperties.MotorDist = [];
    end
    
else
    [cellProperties.Inputs, cellProperties.PSDID] = SynapticPartners(cellID,synapseIden,df);
    
    for i =1:size(cellProperties.PSDID,1)
        cellProperties.PSDsize(i,1) = df.size(df.psd_segid == cellProperties.PSDID(i));
    end
    
    cellProperties.PreSynCoords = PrePartnerCoordinates(cellProperties.PSDID,df);
    cellProperties.PreSynCoordsTransformed = TransformPoints(cellProperties.PreSynCoords,0);
    
    if (isExistReRoot(cellID) ==1 )
        cellProperties.PathLength =  PathLengthToCoordinate(cellProperties.PreSynCoordsTransformed,cellProperties.Tree{1});
    else
        cellProperties.PathLength = [];
    end
    
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
    
    cellProperties.MotorDist = isMotor(cellID,df); % 1X5 (cellID, ABDr, ABDc, ABDir, ABDic)
    
end

%cellProperties.InputsRhombomeres = struct2array(isRhombomere(cellProperties.Inputs));


end
