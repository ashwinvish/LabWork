function [eyeVariables] = eyeVars(eyePosition, time)
% eyeBehaviorCorrelations calculates the Pearson coorelation of the Calcium
% traces to the eyePositon
% baselineCorrectedCalcium , baseline is calculate from
% regressWithBaseline.m
% eyePositon is the trace of the eye positions 
% The behavioral variables are 
% [eyePositionFiltered eyePositionLeft eyePositionRight eyeVelocityLeft eyeVelocityRight] 

% time steps
delt=mean(diff(time));
% filter eye position
eyePositionFiltered=medfilt1(eyePosition,51);
eyePositionFiltered=eyePositionFiltered';
% calculate eye velocity
eyeVelocity=[0;diff(eyePositionFiltered)/delt];
eyeVelocityLeft=eyeVelocity;
eyeVelocityRight=-eyeVelocity;
eyeVelocityLeft(eyeVelocityLeft<0)=0;
eyeVelocityRight(eyeVelocityRight<0)=0;
% calculate left and Right eye positions
eyePositionLeft=eyePositionFiltered;
eyePositionRight=-eyePositionFiltered;
eyePositionLeft(eyePositionLeft<0)=0;
eyePositionRight(eyePositionRight<0)=0;


eyeVariables=[eyePositionFiltered eyePositionLeft eyePositionRight eyeVelocityLeft eyeVelocityRight];
eyeVariables = corr(baselineCorrectdCalcium, eyeVariables,'Type','pearson');

end


