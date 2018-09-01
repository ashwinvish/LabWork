function [calciumTraceBaselineCorrected,eyeFit,baseline,rsq,cor1,cor2,beta]=regressWithBaseline(rawCalcuimTrace,eyePos,time,regressors);
% regressWithBaseline calculated the regression coeffecient as listed in
% Miri et al. and add a regularizer as listed in Daie et. al. It is based
% on baseline_regress.
% rawCalciumTrace is the interpolated time corrected calcium trace
% eyePos is the interpolated eye positions
% time is the duration of the imaging run
% regressors are 'position', 'velocity' or 'both'
% % baseline is the baseline that is subtracted from the rawCalcuimTrace to
% give calciumTraceBaselineCorrected
% eyeFit is the fitted activity traces
% rqs is the root mean square error between calciumTraceBaselineCorrected
% and eyeFit

timeDiff=mean(diff(time));
eyePos=medfilt1(eyePos,51); % smooth eye position signals
eyePos=eyePos';

eyeVelocity=[0;diff(eyePos)/timeDiff]; % calculate eye velocity
eyeVelocityLeft=eyeVelocity;
eyeVelocityRight=-eyeVelocity;
eyeVelocityLeft(eyeVelocityLeft<0)=0;
eyeVelocityRight(eyeVelocityRight<0)=0;

eyePosLeft=eyePos;
eyePosRight=-eyePos;
eyePosLeft(eyePosLeft<0)=0;
eyePosRight(eyePosRight<0)=0;

x=linspace(-1,1,length(eyePos))';
K=ones(size(eyePos));


% compse the linear regression matrix with regressors
if isequal(regressors,'velocity')
    xmat=[eyeVelocity eyeVelocityLeft eyeVelocityRight];
elseif isequal(regressors,'position')
    xmat=[eyePos eyePosLeft eyePosRight];
else
    xmat=[eyePos eyePosLeft eyePosRight eyeVelocityLeft eyeVelocityRight];
end

xmat=xmat-repmat(xmat(1,:),size(xmat,1),1);

for i=1:size(xmat,2);
    a=conv(xmat(:,i),exp(-time/1.9)); % convolve with the CIRF kernel. Time constant is 1.9 sec
    xmat(:,i)=a(1:length(time));
end
stp=size(xmat,2)+1;
xmat=[xmat K x x.^2 x.^3]; % add term to regularize

for i=1:size(rawCalcuimTrace,2);
    beta=pinv(xmat)*rawCalcuimTrace(:,i);
    fit(:,i)=xmat*beta;
end

baseline=xmat(:,stp:end)*beta(stp:end);
eyeFit=xmat(:,1:stp-1)*beta(1:stp-1);
calciumTraceBaselineCorrected=rawCalcuimTrace-baseline;
beta=beta(stp:end);
rsq=1-sum((calciumTraceBaselineCorrected-eyeFit).^2)/sum((calciumTraceBaselineCorrected-mean(calciumTraceBaselineCorrected)).^2);
k=sort(calciumTraceBaselineCorrected);
mn=k(floor(length(k)*.2)); % find the bottom 20% of intensities
a = find(calciumTraceBaselineCorrected<mn);
mn = mean(calciumTraceBaselineCorrected(a));
calciumTraceBaselineCorrected = calciumTraceBaselineCorrected - mn;
cor1=corr(calciumTraceBaselineCorrected,eyeFit,'type','spearman');

if isequal(regressors,'velocity')
    cor2 = corr(calciumTraceBaselineCorrected,xmat(:,1:3),'type','spearman');
elseif isequal(regressors,'position')
    cor2 = corr(calciumTraceBaselineCorrected,xmat(:,1:3),'type','spearman');
else
    cor2 = corr(calciumTraceBaselineCorrected,xmat(:,1:5),'type','spearman');
end

%k=sort(eyfit);mn=k(floor(length(k)*.2));a=find(eyfit<mn);
%mn=mean(cl(a));
%cl=cl-mn;