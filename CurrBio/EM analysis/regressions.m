clear
addpath(genpath('/Users/admin/Documents/spikes'));
addpath(genpath('/Users/admin/Documents/brick'));
addpath(genpath('/Users/admin/Documents/Scripts'));

% load data
load('101112 _files1_4.mat');   % var name is SPT
load('101112_NewROIs.mat');     % var name is NewROIs

% construct new structure

SPTnew = SPT;
for kk = 1:3
    planeNo = kk+1;
    Corig = mat2cell(SPT(planeNo).intensity(:,:),size(SPT(planeNo).intensity(:,:),1), ones(1,size(SPT(planeNo).intensity(:,:),2)));
    Cnew  = mat2cell(NewROIs(planeNo).intensity(:,:),size(NewROIs(planeNo).intensity(:,:),1), ones(1,size(NewROIs(planeNo).intensity(:,:),2)));
    Cmerge = horzcat(Corig,Cnew);
    ImagedTimes = [SPT(planeNo).time,NewROIs(planeNo).time];
    SPTnew(planeNo).numROIs = size(Cmerge,2);
    SPTnew(planeNo).centroid = [SPT(planeNo).centroid, NewROIs(planeNo).centroids];
    SPTnew(planeNo).time = [SPT(planeNo).time, NewROIs(planeNo).time];
    SPTnew(planeNo).intensity = [SPT(planeNo).intensity , NewROIs(planeNo).intensity];
    
    % calulate correlations, use the ipsi eye for corelations; SPT(2).theta(:,3)
    for i = 1:size(Cmerge,2)
        SPTnew(planeNo).cor(1:5,i) = corr(SPTnew(planeNo).intensity(1:end-2,i),SPT(2).eyevars,'Type','Pearson')';
    end
end

% % calculate imaged time
% a=[SPT(planeNo).time(end,:) SPT(planeNo).theta(end,1)];
% a(a<2)=[];
% tstr=max([SPT(planeNo).time(1,:) SPT(planeNo).theta(1,1)]);
% tend=min(a);
% delt=.05; % interpolating steps
% t=tstr:delt:tend;
% 
% 
% % interpolated fluorescence
% 
% for kk = 1:size(Cnorm,2) 
%    Cinterp{kk} = interp1(ImagedTimes(:,kk),sp(:,kk),t,'linear');
% end
% 
% for i = 1:size(Cinterp,2)
%     [CintepBaselineCorrected{i},eyeFits{i},baseline{i},Rsqs(i),Corr(i)]=baseline_regress(Cinterp{i}',e,t,e*0);
% end

%% Compute saccade triggered averages
FLUOR(1).ROI=SPTnew;
FLUOR=stimSummaryEM(FLUOR);
%%
delt=.05;
STA=FLUOR(1).STA;
SPT=FLUOR(1).ROI;
sta=[STA.staf];
err=[STA.err];
er=(sta./err);
er=mean(er);
rsqs=[SPT.rsqs];
for i=1:size(sta,2);
    confplot(sta(:,i),err(:,i));
    title(num2str(er(i)));
    pause;
end
staE=[STA.staE];
stc=diag(corr(sta(1:80,:),staE(1:80,:)))';
load('calibrate1');
cor=max(abs([SPT.cor]));
% a=find(er>3 & cor>.2);
ind=size(SPT(1).intensity,2);
% for i=2:length(SPT);
%     areaI=[areaI SPT(i).areaI+ind];
%     b=size(SPT(i).intensity,2);
%     ind=ind+b;
% end
x=[SPT.centroid];
z=[SPT.z];
int=[SPT.intCorrect];

a=find(stc>.6);
a=intersect(a,find(rsqs>.25));
a=intersect(a,find(cor>.2));
a=sort([a size(SPT(1).intensity,2)+[17 9]]);
% a=intersect(a,areaI);
a=intersect(a,find(z~=0));
sta=sta(22:end,a);
x=x(:,a);x/512*calibrate1(3);
x=[x;z(a)];
n=length(a);
int=int(:,a);
for i=1:n;
    for j=1:n;
        d(i,j)=sqrt(sum((x(:,i)-x(:,j)).^2));
    end
end

load('modes_stafs');
fs=[fs sta];
t=0:delt:delt*(size(fs,1)-1);
t=t-2;
r=svd_rates(fs,t,1.9,[3,3],-100);
r=r(:,end-(n-1):end);
firing=r;
fluorescence=sta(1:end-1,:);
T=t(1:end-1);
r=r(40:end,:);
t=0:.05:.05*(size(r,1)-1);
pwCor=corr(int);
o=find(eye(n,n)==0);
% [Xs,Ys,stds,Cs,Ps]=mean_bin_plot(d(o),pwCor(o),5,1,1);
staE=[STA.staE];eyes=mean(staE,2);