clear;
delt = 0.05;
load('modes_stafs');
load('r_original.mat');
load('STA_raw.mat');
load('MelanieDBXCells.mat');

A = STA_raw;
A = [DBX_vglut/100, A];
A = A(22:end,:);
%A = [fs,A];


t=0:delt:delt*(size(A,1)-1);
t=t-2;
r=svd_rates(A,t,1.9,[3,3],-100);
%r=r(:,end-(n-1):end);
firing=r;
t=0:.05:.05*(size(r,1)-1);

eva = evalclusters(r','kmeans','silhouette','Distance','correlation','KList',[1:6]);
plot(eva);

figure;

shadedErrorBar(t,mean(r(:,eva.OptimalY==1),2), std(r(:,eva.OptimalY==1),[],2));
hold on;
shadedErrorBar(t,mean(r(:,eva.OptimalY==2),2), std(r(:,eva.OptimalY==2),[],2));


%% deconvolution
tau_cirf = 1.9;
dlen = floor(5/delt);       % length of kernel (in index units)
Tnew = linspace(0,dlen*delt,dlen);
ker = exp(-(Tnew-Tnew(1))/tau_cirf);
zpad = zeros(dlen-1,1);     % zero-padding

for i = 1:size(r,2)
    deconvA(:,i) = (tau_cirf/delt)*deconv([A(:,i);zpad],ker);
    medfilt1(deconvA(:,i));
end
