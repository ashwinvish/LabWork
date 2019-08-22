for i = 1:numel(ALX.cellIDs)
    index = find(ismember([Int.cellID]',ALX.cellIDs(i)));
    if ~isempty(index)
        ALX.SaccadicNumbers(i) = length(Int(index).Saccadic);
        ALX.VestibularNumbers(i) = length(Int(index).Vestibular);
        ALX.IntegratorNumbers(i) = length(Int(index).Integrator);
        ALX.ContraNumbers(i) = length(Int(index).Contra);
        ALX.RestNumbers(i) = length(Int(index).EverythingElse);
        if isempty(Int(index).Origin)
            ALX.Origin(i,:) = [NaN, NaN, NaN];
        else
            ALX.Origin(i,:) = Int(index).Origin;
        end
    end
end