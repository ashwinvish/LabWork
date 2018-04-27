clear
addpath(genpath('/Users/admin/Documents/spikes'));
addpath(genpath('/Users/admin/Documents/brick'));
addpath(genpath('/Users/admin/Documents/Scripts'));

%% Load data
load('101112 _files1_4.mat');   % var name is SPT
load('101112_NewROIs.mat');     % var name is NewROIs
i = 2;
Corig = mat2cell(SPT(i).intensity(:,:),size(SPT(i).intensity(:,:),1), ones(1,size(SPT(i).intensity(:,:),2)));
Cnew  = mat2cell(NewROIs(i).intensity(:,:),size(NewROIs(i).intensity(:,:),1), ones(1,size(NewROIs(i).intensity(:,:),2)));
Cmerge = horzcat(Corig,Cnew);
Cmean =cellfun(@mean,Cmerge);
Cmean = num2cell(Cmean);
Cnorm = cellfun(@rdivide,Cmerge,Cmean,'UniformOutput',0); % divide by the mean

a=[SPT(i).time(end,:) SPT(i).theta(end,1)];
a(a<2)=[];
strt=max([SPT(i).time(1,:) SPT(i).theta(1,1)]);
tend=min(a);
delt=.05; % interpolating steps
t=strt:delt:tend;

% interpolate eye position
E=interp1(SPT(i).theta(:,1),SPT(i).theta(:,3),t,'linear');
k=sort(E);
l=length(k);
bot=k(floor(l*.05));
top=k(floor(l*.95));
E=E-((top-bot)/2+bot);
e=(E)*30/(top-bot);

ImagedTimes = [SPT(i).time,NewROIs(i).time];
sp=cell2mat(Cmerge);

for kk = 1:size(Cnorm,2) 
   Cinterp{kk} = interp1(ImagedTimes(:,kk),sp(:,kk),t,'linear');
end

% get eyemovement

Eposmean = e;
Eposmean = 0.05*Eposmean; % scale Eposmean

Epos = cell(1,40);
Epos(:) = {e};

% consolidate all cells



% display image
figure(1234)
imagesc(SPT(i).mnImage)
hold on;
scatter(SPT(i).centroid(1,:), SPT(i).centroid(2,:),'go','filled');
text(10+SPT(i).centroid(1,:), 10+SPT(i).centroid(2,:),num2cell(1:size(SPT(i).centroid(2,:),2)),'color','g');
scatter(NewROIs(i).centroids(1,:), NewROIs(i).centroids(2,:),'ro','filled');
text(10+NewROIs(i).centroids(1,:), 10+NewROIs(i).centroids(2,:),...
num2cell(size(SPT(i).centroid(2,:),2)+1:1:size(SPT(i).centroid(2,:),2)+size(NewROIs(i).centroids(2,:),2)),'color','r');
axis image;
axis off;
set(gca,'XDir','reverse');
colormap gray;


%% run ML spike 

dt = 0.9766; % 
% parameters
% (get the default parameters)
par = tps_mlspikes('par');
% (indicate the frame duration of the data)
par.dt = dt;
% (set physiological parameters)
par.a = 0.091; % DF/F for one spike
par.tau = 0.6539; % 1.936 decay time constant (second)
par.saturation = 0.1; % OGB dye saturation
% (set noise parameters)
par.finetune.sigma = .0379; % a priori level of noise (if par.finetune.sigma
                          % is left empty, MLspike has a low-level routine 
                          % to try estimating it from the data
                          %[0.0334,0.0404,0.0389,0.0389]
par.drift.parameter = .01; % if par.drift parameter is not set, the 
                           % algorithm assumes that the baseline remains
                           % flat; it is also possible to tell the
                           % algorithm the value of the baseline by setting
                           % par.F0
% (do not display graph summary)
par.dographsummary = false;

% spike estimation
[spikest fit drift] = spk_est(Cinterp,par);

% display
figure(1)
spk_display(dt,{spikest},{Cinterp fit drift Epos},'title',num2cell(1:size(Cmerge,2)),'paper');
set(1,'numbertitle','off','name','MLspike alone');
%%
figure(2)
eyeFits(fit,e,Cinterp,t,[]);
%% Auto-calibration
% Now we want to estimate the original spikes from the calcium signals
% rather than using fixed parameter.
% First we run the autocalibration algorithm to attempt estimating, from
% the calcium signals directly, the values of parameters A, tau and sigma.

% parameters
% (get default parameters and set frame duration)
pax = spk_autocalibration('par');
pax.dt = 0.05;
% (set limits for A and tau)
amin = 0.05;
amax = 0.1;
pax.amin = amin;
pax.amax = amax;
a = amin * exp(rand(1)*log(amax/amin)); % log-uniform distribution 
taumin =  0.4;
taumax = 2;
pax.taumin = taumin;
pax.taumax = taumax;
tau = taumin * exp(rand(1)*log(taumax/taumin));

% (set saturation parameter)
pax.saturation = 0.1;
% (give real values to the algorithm for display purposes - they obviously won't be used in the estimation!)
%pax.realspikes = spikes;

pax.reala = a;
pax.realtau = tau;
% (when running MLspike from spk_autocalibratio, do not display graph summary)
pax.mlspikepar.dographsummary = false;

% perform auto-calibration
[tauest aest sigmaest] = spk_autocalibration(Cinterp(1:10),pax)