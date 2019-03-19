function [popMotorDiff] = MotorDiff(cellID,df);

for i = 1:numel(cellID)
    MotorDistribution(i,:) = isMotor(cellID(i),df);
end

popMotorDiff = (MotorDistribution(:,2)+MotorDistribution(:,3)) - ...
        (MotorDistribution(:,4)+MotorDistribution(:,5));
    
end