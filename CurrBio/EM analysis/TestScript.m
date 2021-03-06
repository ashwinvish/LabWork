clear

if ismac
    addpath(genpath('/Users/admin/Documents/spikes'));
    addpath(genpath('/Users/admin/Documents/brick'));
    addpath(genpath('/Users/admin/Documents/Scripts'));
else
    addpath(genpath('/usr/people/ashwinv/seungmount/research/Ashwin/spikes'));
    addpath(genpath('/usr/people/ashwinv/seungmount/research/Ashwin/brick'));
    addpath(genpath('/usr/people/ashwinv/seungmount/research/Ashwin/Scripts'));
end

%% Load data
load('101112 _files1_4.mat');   % var name is SPT
load('101112_NewROIs.mat');     % var name is NewROIs
imagingPlane = 4;
eyePlane = [2,3,2,3]; % the eye that contained the best SNR was considered
Corig = mat2cell(SPT(imagingPlane).intensity(:,:),size(SPT(imagingPlane).intensity(:,:),1), ones(1,size(SPT(imagingPlane).intensity(:,:),2)));
Cnew  = mat2cell(NewROIs(imagingPlane).intensity(:,:),size(NewROIs(imagingPlane).intensity(:,:),1), ones(1,size(NewROIs(imagingPlane).intensity(:,:),2)));
Cmerge = horzcat(Corig,Cnew);
Cmean =cellfun(@mean,Cmerge);
Cmean = num2cell(Cmean);
Cnorm = cellfun(@rdivide,Cmerge,Cmean,'UniformOutput',0); % divide by the mean

a=[SPT(imagingPlane).time(end,:) SPT(imagingPlane).theta(end,1)];
a(a<2)=[];
strt=max([SPT(imagingPlane).time(1,:) SPT(imagingPlane).theta(1,1)]);
tend=min(a);
delt=.05; % interpolating steps
interpolatedTime=strt:delt:tend;

% interpolate eye position
E=interp1(SPT(imagingPlane).theta(:,1),SPT(imagingPlane).theta(:,eyePlane(imagingPlane)),interpolatedTime,'linear');
k=sort(E);
l=length(k);
bot=k(floor(l*.05));
top=k(floor(l*.95));
E=E-((top-bot)/2+bot);
eyePositon=(E)*30/(top-bot);

ImagedTimes = [SPT(imagingPlane).time,NewROIs(imagingPlane).time];
rawCalciumTraces=cell2mat(Cmerge);

for kk = 1:size(Cnorm,2) 
   Cinterp{kk} = interp1(ImagedTimes(:,kk),rawCalciumTraces(:,kk),interpolatedTime,'linear');
end

% get eyemovement

Eposmean = eyePositon;
Eposmean = 0.05*Eposmean; % scale Eposmean

Epos = cell(1,40);
Epos(:) = {eyePositon};

% consolidate all cells



% % display image
% figure(1234)
% imagesc(SPT(imagingPlane).mnImage)
% hold on;
% scatter(SPT(imagingPlane).centroid(1,:), SPT(imagingPlane).centroid(2,:),'go','filled');
% text(10+SPT(imagingPlane).centroid(1,:), 10+SPT(imagingPlane).centroid(2,:),num2cell(1:size(SPT(imagingPlane).centroid(2,:),2)),'color','g');
% scatter(NewROIs(imagingPlane).centroids(1,:), NewROIs(imagingPlane).centroids(2,:),'ro','filled');
% text(10+NewROIs(imagingPlane).centroids(1,:), 10+NewROIs(imagingPlane).centroids(2,:),...
% num2cell(size(SPT(imagingPlane).centroid(2,:),2)+1:1:size(SPT(imagingPlane).centroid(2,:),2)+size(NewROIs(imagingPlane).centroids(2,:),2)),'color','r');
% axis image;
% axis off;
% set(gca,'XDir','reverse');
% colormap gray;

centeroids = [SPT(imagingPlane).centroid,NewROIs(imagingPlane).centroids];

figure(1234)
imagesc(SPT(imagingPlane).mnImage')
hold on;
scatter(SPT(imagingPlane).centroid(2,:), SPT(imagingPlane).centroid(1,:),'go','filled');
text(10+SPT(imagingPlane).centroid(2,:), 10+SPT(imagingPlane).centroid(1,:),num2cell(1:size(SPT(imagingPlane).centroid(1,:),2)),'color','g');
scatter(NewROIs(imagingPlane).centroids(2,:), NewROIs(imagingPlane).centroids(1,:),'ro','filled');
text(10+NewROIs(imagingPlane).centroids(2,:), 10+NewROIs(imagingPlane).centroids(1,:),...
num2cell(size(SPT(imagingPlane).centroid(2,:),2)+1:1:size(SPT(imagingPlane).centroid(1,:),2)+size(NewROIs(imagingPlane).centroids(1,:),2)),'color','r');
axis image;
axis off;
%set(gca,'XDir','reverse');
colormap gray;
title(imagingPlane);



%% Old coorelations of activity with eyevars. This baseline correction was used for the CB paper
rawCalciumTracesBaseline = zeros(size(rawCalciumTraces));
x=linspace(-1,1,size(rawCalciumTraces,1))';
xmat=[ones(size(x)) x];
clear x;
for i = 1: size(rawCalciumTraces,2)
    baseline = xmat*pinv(xmat)*rawCalciumTraces(:,i);
    rawCalciumTracesBaseline(:,i)=rawCalciumTraces(:,i)-baseline; %subtract the baseline
end
NewROIs(imagingPlane).intCorrected = rawCalciumTracesBaseline(:,size(Corig,2)+1:end);
NewROIs(imagingPlane).cor = corr(rawCalciumTracesBaseline(1:end-2,size(Corig,2)+1:end),SPT(2).eyevars);

%% New coorelations of activity with eyevars based on supplementary of Daie et.al 2015

for i = 1:size(Cinterp,2)
    %since eyepos recordings have higher temporal resolution, find the
    %appropriate eyepos for given cell
    for jj = 1:size(ImagedTimes,1)
        [~,b2(jj)] = min(abs(ImagedTimes(jj,i)-SPT(imagingPlane).theta(:,1)));
    end
    eyePos = SPT(imagingPlane).theta(b2,eyePlane(imagingPlane)); % get eye positions
    k=sort(eyePos); 
    l=length(k);
    bot=k(floor(l*.05));
    top=k(floor(l*.95));
    eyePos= eyePos-((top-bot)/2+bot); % clean eye positions
    timeDiff=mean(diff(ImagedTimes(:,i))); % compute eye vars
    eyePos=eyePos; 
    eyeVelocity=[0;diff(eyePos)/timeDiff]; % calculate eye velocity
    eyeVelocityLeft=eyeVelocity;
    eyeVelocityRight=-eyeVelocity;
    eyePosRight(eyeVelocityLeft<0)=0;
    eyeVelocityRight(eyeVelocityRight<0)=0;
    
    eyePosLeft=eyePos;
    eyePosRight=-eyePos;
    eyePosLeft(eyePosLeft<0)=0;
    eyePosRight(eyePosRight<0)=0;
    
    eyeVars{i} = [eyePos,eyePosLeft,eyePosRight,eyeVelocityLeft,eyeVelocityRight];
    [CNewBaselineCorrected{i},~,~,RsqsRawAndFit]=regressWithBaseline(rawCalciumTraces(:,i),eyePos',ImagedTimes(:,i)','both');
    eyeVarsCorrNonInterpolated(i,:) = corr(CNewBaselineCorrected{i},eyeVars{i});
    clear b;
end

%% Make plots based ion raw, non interpolated signals

figure;
eyeVarsMat = cell2mat(eyeVars);
meanEyeVars = mean(eyeVarsMat(:,1:5:end),2);
dff = cell2mat(CNewBaselineCorrected);
%dff = (dff-mean(dff))./mean(dff);
m1 = [min(dff(:)), max(dff(:))];
map = colorcet('L8');
%map = parula(256);
%map= cbrewer('seq','PuBu',128);
%map = NLmap(map,m1(2),m1(1),m1(1),0.05);
plottingTime = mean(ImagedTimes,2);

subplot(4,1,1)
plot(plottingTime, meanEyeVars);
xlim([0, max(plottingTime)]);
xticks(plottingTime(50:50:300));
ylabel('eyePositon');
box off;
title(imagingPlane);

subplot(4,1,2)
imagesc(dff');
xticklabels(plottingTime(50:50:300));
ylabel('cellnumber');
title('no sorting');
%xlabel('time');
box off;
colormap(map)

% correlated to eyepositions
maxEyeCorrNonInterpolated_P = max(eyeVarsCorrNonInterpolated(:,1:3),[],2);
[a1,b1] = sort(maxEyeCorrNonInterpolated_P,'descend');

subplot(4,1,3);
imagesc(dff(:,b1)');
xticklabels(plottingTime(50:50:300));
%yyaxis left
ylabel('cellnumber');
yticks(2:2:size(b1,1))
yticklabels(b1(2:2:size(b1,1)));
% yyaxis right
% yticks(1:4:size(b1,1))
% yticklabels(b1(1:4:size(b1,1)));
%ytickangle(45);
title('sorted by coorelation to eyePositon');
%xlabel('time');
box off;
colormap(map);

maxEyeCorrNonInterpolated_V = max(eyeVarsCorrNonInterpolated(:,4:5),[],2);
[a2,b2] = sort(maxEyeCorrNonInterpolated_V,'descend');

subplot(4,1,4);
imagesc(dff(:,b2)');
xticklabels(plottingTime(50:50:300));
ylabel('cellnumber');
yticks(2:2:size(b2,1))
yticklabels(b2(2:2:size(b2,1)));
%ytickangle(45);
title('sorted by coorelation to eye velocity');
%xlabel('time');
box off;
colormap(map)

figure

subplot(2,2,1)
imagesc(SPT(imagingPlane).mnImage')
hold on;
scatter(SPT(imagingPlane).centroid(2,:), SPT(imagingPlane).centroid(1,:),'go','filled');
text(10+SPT(imagingPlane).centroid(2,:), 10+SPT(imagingPlane).centroid(1,:),num2cell(1:size(SPT(imagingPlane).centroid(1,:),2)),'color','g');
scatter(NewROIs(imagingPlane).centroids(2,:), NewROIs(imagingPlane).centroids(1,:),'ro','filled');
text(10+NewROIs(imagingPlane).centroids(2,:), 10+NewROIs(imagingPlane).centroids(1,:),...
num2cell(size(SPT(imagingPlane).centroid(2,:),2)+1:1:size(SPT(imagingPlane).centroid(1,:),2)+size(NewROIs(imagingPlane).centroids(1,:),2)),'color','r');
axis image;
axis off;
%set(gca,'XDir','reverse');
colormap gray;
title(imagingPlane);


subplot(2,2,3)
%colorMap = cbrewer('seq','RdPu',40)
%ax1 = axes;
scatter(centeroids(2,b1), centeroids(1,b1),50,maxEyeCorrNonInterpolated_P(b1),'o','filled');
caxis([min(maxEyeCorrNonInterpolated_P),max(maxEyeCorrNonInterpolated_P)]);
colorbar;
set(gca,'YDir','reverse','color','none','XLim',[0.5,512.5],'YLim',[0.5,512.5]);
xlabel('M <------> L');
ylabel('C <------> R');
title('correlated to eye position');
axis on;
box on;
daspect([1,1,1]);

subplot(2,2,4)
%colorMap = parula(size(b2,1));
scatter(centeroids(2,b2), centeroids(1,b2),50,maxEyeCorrNonInterpolated_V(b2),'o','filled');
caxis([min(maxEyeCorrNonInterpolated_V),max(maxEyeCorrNonInterpolated_V)]);
colorbar;
set(gca,'YDir','reverse','color','none','XLim',[0.5,512.5],'YLim',[0.5,512.5]);
axis off;
title('correlated to eye velocity');
xlabel('M <------> L');
ylabel('C <------> R');
axis on;
box on
daspect([1,1,1]);


figure;
scatter(maxEyeCorrNonInterpolated_P(b1),maxEyeCorrNonInterpolated_V(b2),'ko','filled');
hold on;
line([-0.05,6],[-0.05,6],'color','k','linestyle','--');
%line([-0.05,0.2],[0.2,0.2], 'color','k');
%line([0.2,0.2],[-0.05,0.2], 'color','k');
set(gca, 'XLim',[-0.05,0.6],'YLim',[-0.05,0.6]);
xlabel('Correlation to eye position');
ylabel('Correlation to eye velocity');
axis square;



%% run baseline_regress to obtain the fits
% for this process the raw calcium traces need to be converted to df/f
% for i = 1:size(Cinterp,2)
%     Cinterpdff{i} = (Cinterp{i}-mean(Cinterp{i}))/mean(Cinterp{i});
% end

for i = 1:size(Cinterp,2)
    [CinterpBaselineCorrected_V{i},eyeFit_V{i},baseline,Rqs_V(i),baselineToFitCor_V(i),eyeVarCorInterpolated_V(:,i),beta]=regressWithBaseline(Cinterp{i}',eyePositon,interpolatedTime,'velocity');
    [CinterpBaselineCorrected_P{i},eyeFit_P{i},baseline,Rqs_P(i),baselineToFitCor_P(i),eyeVarCorInterpolated_P(:,i),beta]=regressWithBaseline(Cinterp{i}',eyePositon,interpolatedTime,'position');
    [CinterpBaselineCorrected_B{i},eyeFit_B{i},baseline,Rqs_B(i),baselineToFitCor_B(i),eyeVarCorInterpolated_B(:,i),beta]=regressWithBaseline(Cinterp{i}',eyePositon,interpolatedTime,'both');
end

eyeThr=[3 140];
fixationnLength=7;
planeNumber = 2;
STA_V = spontSTtrials(planeNumber,interpolatedTime,cell2mat(CinterpBaselineCorrected_V),eyePositon,eyeThr,fixationnLength);
STA_P = spontSTtrials(planeNumber,interpolatedTime,cell2mat(CinterpBaselineCorrected_P),eyePositon,eyeThr,fixationnLength);
STA_B = spontSTtrials(planeNumber,interpolatedTime,cell2mat(CinterpBaselineCorrected_B),eyePositon,eyeThr,fixationnLength);

%% Plot saccade triggered params based on interpolated signal
figure;
% coorelation to position
for i = 1:size(Cinterp,2)
    subplot(8,8,i)
    yyaxis left; 
    shadedErrorBar(STA_P.T, STA_P.staf(:,i),STA_P.err(:,i),'b');
    yyaxis right; plot(STA_P.T, STA_P.staE(:,i));
    xlim([-2,7]);
    %legend({'Vel','Pos','both'});
    titleName = sprintf('%2d,%1.3f,%1.3f,%1.3f',i,maxEyeCorrNonInterpolated_P(i),...
        maxEyeCorrNonInterpolated_V(i), STA_P.cirfTau(i));
    title(titleName,'FontSize',6);
    box off;
end


% figure; % coorelation to velocity
% for i = 1:size(Cinterp,2)
%     subplot(8,8,i)
%     yyaxis left; 
%     shadedErrorBar(STA_V.T, STA_V.staf(:,i),STA_V.err(:,i),'r');
%     hold on;
%     yyaxis right; plot(STA_V.T, STA_V.staE(:,i));
%     titleName = sprintf('Id:%2d,X_V: %1.3f,T_V:%1.3f',i,STA_V.cor(i),STA_V.cirfTau(i));
%     title(titleName,'FontSize',6);
% end

%% Make plots based on interpolated Ca signals

% eyeVarsMat = cell2mat(eyeVars);
% meanEyeVars = mean(eyeVarsMat(:,1:5:end),2);
% dffIntep = cell2mat(CinterpBaselineCorrected_B);
% dffIntep = cell2mat(eyeFit_B);
% 
% 
% dffIntep = (dffIntep-mean(dffIntep))./mean(dffIntep);
% m1 = [min(dffIntep(:)), max(dffIntep(:))];
% map = colorcet('L8');
% map = parula(256);
% map= cbrewer('seq','PuBu',128);
% map = NLmap(map,m1(2),m1(1),m1(1),60);
% 
% figure;
% subplot(4,1,1)
% plot(interpolatedTime, eyePositon);
% ylabel('eyePositon');
% xlim([1 max(interpolatedTime)]);
% box off;
% title(imagingPlane);
% 
% subplot(4,1,2)
% imagesc(dffIntep');
% ylabel('cellnumber');
% xticks([]);
% title('no sorting');
% xlabel('time');
% box off;
% colormap(map)
% 
% correlated to eyepositions
% maxEyeCorrNonInterpolated_P = max(eyeVarCorInterpolated_B(1:3,:));
% [a3,b3] = sort(maxEyeCorrNonInterpolated_P,'descend');
% 
% subplot(4,1,3);
% imagesc(dffIntep(:,b3)');
% xticks([]);
% ylabel('cellnumber');
% yticks(2:4:size(b3,1))
% yticklabels(b3(2:4:size(b3,1)));
% ytickangle(45);
% title('sorted by coorelation to eyePositon');
% xlabel('time');
% box off;
% colormap(map)
% 
% maxEyeCorrNonInterpolated_V = max(eyeVarsCorrNonInterpolated(4:5,:));
% [a4,b4] = sort(maxEyeCorrNonInterpolated_V,'descend');
% 
% subplot(4,1,4);
% imagesc(dffIntep(:,b4)');
% ylabel('cellnumber');
% yticks(2:4:size(b4,1))
% yticklabels(b4(2:4:size(b4,1)));
% ytickangle(45);
% title('sorted by coorelation to eye velocity');
% xlabel('time');
% box off;
% colormap(map)
% 
% figure;
% subplot(1,2,1)
% colorMap = parula(size(b3,1));
% ax1 = axes;
% scatter(centeroids(2,b3), centeroids(1,b3),50,colorMap,'o','filled');
% colorbar
% set(gca,'YDir','reverse','color','none','XLim',[0.5,512.5],'YLim',[0.5,512.5]);
% title('correlated to eye position');
% xlabel('M <------> L');
% ylabel('C <------> R');
% axis on;
% box on;
% daspect([1,1,1]);
% 
% subplot(1,2,2)
% colorMap = parula(size(b4,1));
% scatter(centeroids(2,b4), centeroids(1,b4),50,colorMap,'o','filled');
% colorbar;
% set(gca,'YDir','reverse','color','none','XLim',[0.5,512.5],'YLim',[0.5,512.5]);
% axis off;
% title('correlated to eye velocity');
% xlabel('M <------> L');
% ylabel('C <------> R');
% axis on;
% box on
% daspect([1,1,1]);
% 
% figure;
% scatter(maxEyeCorrP(b3),maxEyeCorrV(b4),'filled');
% hold  on;
% text(maxEyeCorrP(b1),maxEyeCorrV(b2),num2str(b1),'color','r');
% text(0.01+maxEyeCorrP(b3),maxEyeCorrV(b4),num2str(b4),'color','b');
% ylim([-0.1,max(maxEyeCorrP(b3))]);
% xlim([-0.1,max(maxEyeCorrP(b3))]);
% line([-0.1,max(maxEyeCorrP(b3))], [-0.1,max(maxEyeCorrP(b3))],'color','k','linestyle','--');
% scatter(maxEyeCorrP(b3(find(maxEyeCorrP(b3)>0.2))),maxEyeCorrV(b3(find(maxEyeCorrP(b3)>0.2))),'filled');
% scatter(maxEyeCorrP(b4(find(maxEyeCorrV(b4)>0.2))), maxEyeCorrV(b4(find(maxEyeCorrV(b4)>0.2))),'filled');
% xlabel('Correlation to eye position');
% ylabel('Correlation to eye velocity');
% axis square;
% title(imagingPlane);
% 

%%

    

% %% plots
% figure(2)
% eyeFits(fitting,eyePositon,Cinterp,interpolatedTime);
% %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % run ML spike 
% 
% dt = 0.9766; % 
% % parameters
% % (get the default parameters)
% par = tps_mlspikes('par');
% % (indicate the frame duration of the data)
% par.dt = dt;
% % (set physiological parameters)
% par.a = 0.091; % DF/F for one spike
% par.tau = 0.6539; % 1.936 decay time constant (second)
% par.saturation = 0.1; % OGB dye saturation
% % (set noise parameters)
% par.finetune.sigma = .0379; % a priori level of noise (if par.finetune.sigma
%                           % is left empty, MLspike has a low-level routine 
%                           % to try estimating it from the data
%                           %[0.0334,0.0404,0.0389,0.0389]
% par.drift.parameter = .01; % if par.drift parameter is not set, the 
%                            % algorithm assumes that the baseline remains
%                            % flat; it is also possible to tell the
%                            % algorithm the value of the baseline by setting
%                            % par.F0
% % (do not display graph summary)
% par.dographsummary = false;
% 
% % spike estimation
% [spikest fitting drift] = spk_est(Cinterp,par);
% 
% % display
% % figure(1)
% % spk_display(dt,{spikest},{Cinterp fitting drift Epos},'title',num2cell(1:size(Cmerge,2)),'paper');
% % set(1,'numbertitle','off','name','MLspike alone');
% 
% %% Auto-calibration
% % Now we want to estimate the original spikes from the calcium signals
% % rather than using fixed parameter.
% % First we run the autocalibration algorithm to attempt estimating, from
% % the calcium signals directly, the values of parameters A, tau and sigma.
% 
% % parameters
% % (get default parameters and set frame duration)
% pax = spk_autocalibration('par');
% pax.dt = 0.05;
% % (set limits for A and tau)
% amin = 0.05;
% amax = 0.1;
% pax.amin = amin;
% pax.amax = amax;
% a = amin * exp(rand(1)*log(amax/amin)); % log-uniform distribution 
% taumin =  0.4;
% taumax = 2;
% pax.taumin = taumin;
% pax.taumax = taumax;
% tau = taumin * exp(rand(1)*log(taumax/taumin));
% 
% % (set saturation parameter)
% pax.saturation = 0.1;
% % (give real values to the algorithm for display purposes - they obviously won't be used in the estimation!)
% %pax.realspikes = spikes;
% 
% pax.reala = a;
% pax.realtau = tau;
% % (when running MLspike from spk_autocalibratio, do not display graph summary)
% pax.mlspikepar.dographsummary = false;
% 
% % perform auto-calibration
% [tauest aest sigmaest] = spk_autocalibration(Cinterp(1:10),pax)