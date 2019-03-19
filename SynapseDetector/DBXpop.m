% DBX population distributions
clear;
addpath(genpath('/Users/ashwin/Documents/'));

colors = cbrewer('qual','Dark2',10);
startup

if ismac
    addpath(genpath('/Users/ashwin/Documents/'));
    df = readtable('/Users/ashwin/Documents/SynapseDetector/11252018.csv');
else
    addpath(genpath('/usr/people/ashwinv/seungmount/research/Ashwin/Scripts'));
    df = readtable('/usr/people/ashwinv/seungmount/research/Ashwin/SynapseDetector/11252018.csv');
end


confirmedDBX = [76199 76200 76182 76183 76185 76186 76189 76191 76188 ];
putativeDBX = [77582 77588 77595 77597 77598 77599 77605 79040 76399 76828 ...
    76829 78435 76826 76874 76289 76542 76832 76836 76838 76877 76883 76884 ...
    78400 78429 76383 76867 77594 76830 76876 78440 78438 78447 79045 78434 ...
    78450 76878 77683 77450 78448 ];

AllDBX = [confirmedDBX, putativeDBX];

for i = 1:numel(AllDBX)
    DBX(i) = InputsByClass(AllDBX(i),df);
end

% 
for i = 1:numel(AllDBX)
    DBXSaccadic(i).PathLength = DBX(i).PathLength(DBX(i).isSaccadic)/max(Pvec_tree(DBX(i).Tree{1}));
    DBXVestibular(i).PathLength = DBX(i).PathLength(DBX(i).isVestibular)/max(Pvec_tree(DBX(i).Tree{1}));
    DBXContra(i).PathLength = DBX(i).PathLength(DBX(i).isContra)/max(Pvec_tree(DBX(i).Tree{1}));
    DBXIntegrator(i).PathLength = DBX(i).PathLength(DBX(i).isIntegrator)/max(Pvec_tree(DBX(i).Tree{1}));
    
    DBXSaccadic(i).PSDsize = DBX(i).PSDsize(DBX(i).isSaccadic);
    DBXVestibular(i).PSDsize = DBX(i).PSDsize(DBX(i).isVestibular);
    DBXContra(i).PSDsize = DBX(i).PSDsize(DBX(i).isContra);
    DBXIntegrator(i).PSDsize = DBX(i).PSDsize(DBX(i).isIntegrator)
end


subplot(1,3,1)
[~,plotOrder] = sortrows(vertcat(DBX.Origin),2);
for i = 1:numel(AllDBX)
    scatter(DBXSaccadic(plotOrder(i)).PathLength,DBX(plotOrder(i)).Origin(2)*ones(size(DBXSaccadic(plotOrder(i)).PathLength,1),1),...
        0.05*(DBXSaccadic(plotOrder(i)).PSDsize),colors(1,:));
    hold on;
    %if ~isnan(mean(DBXSaccadic(plotOrder(i)).PathLength))
    %    plot(mean(DBXSaccadic(plotOrder(i)).PathLength),DBX(plotOrder(i)).Origin(2),'-ro','markerSize',0.01*(median(DBXSaccadic(plotOrder(i)).PSDsize)));
    %end
    scatter(DBXVestibular(plotOrder(i)).PathLength,DBX(plotOrder(i)).Origin(2)*ones(size(DBXVestibular(plotOrder(i)).PathLength,1),1),...
        0.05*(DBXVestibular(plotOrder(i)).PSDsize),colors(2,:));
    scatter(DBXContra(plotOrder(i)).PathLength,DBX(plotOrder(i)).Origin(2)*ones(size(DBXContra(plotOrder(i)).PathLength,1),1),...
        0.05*(DBXContra(plotOrder(i)).PSDsize),colors(3,:));
     scatter(DBXIntegrator(plotOrder(i)).PathLength,DBX(plotOrder(i)).Origin(2)*ones(size(DBXIntegrator(plotOrder(i)).PathLength,1),1),...
        0.05*(DBXIntegrator(plotOrder(i)).PSDsize),colors(4,:));
end
%daspect([1,1,1]);
set(gca,'YDir','reverse');
ylabel('RC positon in \mu');
xlabel('Normalized path length');


% subplot(1,3,2);
% [~,plotOrder] = sortrows(vertcat(DBX.Origin),1);
% for i = 1:numel(AllDBX)
%     scatter(DBXSaccadic(plotOrder(i)).PathLength,DBX(plotOrder(i)).Origin(1)*ones(size(DBXSaccadic(plotOrder(i)).PathLength,1),1),...
%         0.05*(DBXSaccadic(plotOrder(i)).PSDsize),'k');
%     hold on;
%     %if ~isnan(mean(DBXSaccadic(plotOrder(i)).PathLength))
%     %    plot(mean(DBXSaccadic(plotOrder(i)).PathLength),DBX(plotOrder(i)).Origin(2),'-ro','markerSize',0.01*(median(DBXSaccadic(plotOrder(i)).PSDsize)));
%     %end
% end
% %daspect([1,1,1]);
% set(gca,'YDir','reverse');
% ylabel('ML positon in \mu');
% xlabel('Normalized path length');
% 
% subplot(1,3,3)
% [~,plotOrder] = sortrows(vertcat(DBX.Origin),3);
% for i = 1:numel(AllDBX)
%     scatter(DBXSaccadic(plotOrder(i)).PathLength,DBX(plotOrder(i)).Origin(3)*ones(size(DBXSaccadic(plotOrder(i)).PathLength,1),1),...
%         0.05*(DBXSaccadic(plotOrder(i)).PSDsize),'k');
%     hold on;
%     %if ~isnan(mean(DBXSaccadic(plotOrder(i)).PathLength))
%     %    plot(mean(DBXSaccadic(plotOrder(i)).PathLength),DBX(plotOrder(i)).Origin(2),'-ro','markerSize',0.01*(median(DBXSaccadic(plotOrder(i)).PSDsize)));
%     %end
% end
% %daspect([1,1,1]);
% set(gca,'YDir','reverse');
% ylabel('DV positon in \mu');
% xlabel('Normalized path length');

%% Saccadic Axons

UniqueDBXSaccadicAxons = unique(vertcat(DBX.Saccadic));
DBXSaccadicMotorDiff = MotorDiff(UniqueDBXSaccadicAxons,df);
histogram(DBXSaccadicMotorDiff,10);

LeadLikeSaccadicAxons = UniqueDBXSaccadicAxons(DBXSaccadicMotorDiff>0);
LagLikeSaccadicAxons = UniqueDBXSaccadicAxons(DBXSaccadicMotorDiff<0);

%transform_swc_AV(LeadLikeSaccadicAxons,colors(1,:),[],true,false);
%transform_swc_AV(LagLikeSaccadicAxons,colors(2,:),[],false,false);
%transform_swc_AV(AllDBX,colors(3,:),[],false,false);

leadOrder = find(ismember(UniqueDBXSaccadicAxons,LeadLikeSaccadicAxons));
lagOrder = find(ismember(UniqueDBXSaccadicAxons,LagLikeSaccadicAxons));

[DBXSaccadicAxons,SaccadicCounts,SaccadicGraph,SaccadicShortestPaths] = PartnerConnectivity(UniqueDBXSaccadicAxons,UniqueDBXSaccadicAxons,df);

SaccadicCounts = SaccadicCounts.*~eye(size(SaccadicCounts));

for i = 1:size(SaccadicShortestPaths,2)
    ii = SaccadicShortestPaths(i).ij(1);
    jj = SaccadicShortestPaths(i).ij(2);
   distMat(ii,jj) = SaccadicShortestPaths(i).dist;
end
distMat = distMat.*~eye(size(distMat));

%% Integrator Axons

UniqueDBXIntegratorAxons = unique(vertcat(DBX.Integrator));
DBXIntegratorMotorDiff = MotorDiff(UniqueDBXIntegratorAxons,df);
histogram(DBXIntegratorMotorDiff,10);

LeadLikeIntegratorAxons = UniqueDBXIntegratorAxons(DBXIntegratorMotorDiff>0);
LagLikeIntegratorAxons = UniqueDBXIntegratorAxons(DBXIntegratorMotorDiff<0);

[DBXIntegratorAxons,IntegratorCounts,IntegratorGraph,IntegratorShortestPaths] = PartnerConnectivity(UniqueDBXIntegratorAxons,UniqueDBXIntegratorAxons,df);

