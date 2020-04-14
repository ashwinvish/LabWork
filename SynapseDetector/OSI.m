function [ocularSelectivity,motorSum,interSum] = OSI(cellIDs,df)
%format longE
motorDistribution = isMotor(cellIDs,df);
motorSum = motorDistribution(:,2)+motorDistribution(:,3);
motorSum = double(motorSum);
interSum = motorDistribution(:,4)+motorDistribution(:,5);
interSum = double(interSum);
ocularSelectivity = (motorSum-interSum)./(motorSum+interSum);
end

