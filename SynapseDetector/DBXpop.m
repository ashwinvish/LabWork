% DBX population distributions
clear;
addpath(genpath('/Users/ashwin/Documents/'));

colors = cbrewer('qual','Dark2',10);

temp1 = colorcet('CBTL1','N',5); %  linear-tritanopic_krjcw_5-98_c46_n256
temp2 = colorcet('CBL2','N',5);  %  linear-protanopic-deuteranopic_kbw_5-98_c40_n256

leadColor = temp1(3,:); % dark red, CB safe
lagColor = temp2(3,:); % dark blue, CB safe

PartnerColors = colorcet('CBTD1','N',10);
lightRed = PartnerColors(10,:);
lightBlue = PartnerColors(1,:);


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
figure;
%histogram(DBXSaccadicMotorDiff,10);

LeadLikeSaccadicAxons = UniqueDBXSaccadicAxons(DBXSaccadicMotorDiff>0);
LagLikeSaccadicAxons = UniqueDBXSaccadicAxons(DBXSaccadicMotorDiff<0);

leadAxonOrder = find(ismember(UniqueDBXSaccadicAxons,LeadLikeSaccadicAxons));
lagAxonOrder = find(ismember(UniqueDBXSaccadicAxons,LagLikeSaccadicAxons));


LeadLikeDBX = []
LagLikeDBX = [];
for i = 1:numel(AllDBX)
    [a,~] = SynapticPartners(AllDBX(i),1,df);
    a = a(a<1e5);
    leadSum = 0;
    lagSum = 0;
    for jj = 1:length(a)
        leadSum = leadSum + sum(a(jj) == LeadLikeSaccadicAxons);
        lagSum = lagSum + sum(a(jj) == LagLikeSaccadicAxons);
    end
    if (leadSum-lagSum) > 0 
        LeadLikeDBX = [LeadLikeDBX; AllDBX(i)];
    else
        LagLikeDBX = [LagLikeDBX; AllDBX(i)];
    end
end


for i = 1:numel(UniqueDBXSaccadicAxons)
    [a,~] = SynapticPartners(UniqueDBXSaccadicAxons(i),2,df);
    a = a(a<1e5);
    leadSynapses = 0;
    lagSynapses = 0;
    for jj = 1:length(a)
    leadSynapses = leadSynapses + sum(a(jj) == LeadLikeDBX);
    lagSynapses = lagSynapses + sum(a(jj) == LagLikeDBX);
    end
    DBXSynapticDiff(i) = leadSynapses - lagSynapses;
end

temp = isPostSynapseMotor(UniqueDBXSaccadicAxons,df);
DBXSaccadicMotorNeuronDiff = sum(temp(:,1:2),2) - sum(temp(:,3:4),2);

MarkerColorOrder(leadAxonOrder,:) = repmat(leadColor,size(leadAxonOrder,1),1);
MarkerColorOrder(lagAxonOrder,:) = repmat(lagColor,size(lagAxonOrder,1),1);

subplot(2,3,3)
scatter(DBXSynapticDiff,DBXSaccadicMotorDiff,20,MarkerColorOrder,'filled');
line([0,0],[-100,100],'linestyle',':','color','k');
line([-10,10],[0,0],'linestyle',':','color','k');
ylabel('#(A - Ai) synapses');
xlabel('difference in synapses')
box off;
axis square;

subplot(2,3,6)

scatter(DBXSaccadicMotorDiff,DBXSaccadicMotorNeuronDiff,20,MarkerColorOrder,'filled');
line([0,0],[-20,20],'linestyle',':','color','k');
line([-200,200],[0,0],'linestyle',':','color','k');
ylabel('#(A - Ai) neurons');
xlabel('#(A - Ai) synapses')
box off;
axis square;

subplot(2,3,[1,4])
transform_swc_AV(LeadLikeDBX,leadColor,[],true,false);
%subplot(2,3,[2,5])
transform_swc_AV(LagLikeDBX,lagColor,[],true,false);

subplot(2,3,[2,5])
LeadDBXOrder = find(ismember(AllDBX,LeadLikeDBX));
LagDBXOrder = find(ismember(AllDBX,LagLikeDBX));
LeadDBXOrigins = vertcat(DBX(LeadDBXOrder).Origin);
LagDBXOrigins = vertcat(DBX(LagDBXOrder).Origin);
histogram(LeadDBXOrigins(:,2),30,'EdgeColor',leadColor,'Orientation','horizontal','DisplayStyle','stairs','LineWidth',2);
hold on;
histogram(LagDBXOrigins(:,2),30,'EdgeColor',lagColor,'Orientation','horizontal','DisplayStyle','stairs','LineWidth',2);
set(gca,'YDir','reverse','YLim',[0,1280]);
daspect([1,30,1]);

%%
[DBXSaccadicAxons,SaccadicCounts,SaccadicGraph,SaccadicShortestPaths] = PartnerConnectivity(UniqueDBXSaccadicAxons,UniqueDBXSaccadicAxons,'Saccadic',df);

SaccadicCounts = SaccadicCounts.*~eye(size(SaccadicCounts));

figure;

subplot(2,3,[1,4])
transform_swc_AV(LeadLikeSaccadicAxons,lightRed,[],true,false);
subplot(2,3,[2,5])
transform_swc_AV(LagLikeSaccadicAxons,lightBlue,[],true,false);

subplot(2,3,3)
imagesc(SaccadicCounts([leadAxonOrder;lagAxonOrder],[leadAxonOrder;lagAxonOrder]))
minMax = [min(SaccadicCounts(:)), max(SaccadicCounts(:))];
colorcet('L1','N',minMax(2)-minMax(1),'reverse',1);
line([0,size(SaccadicCounts,2)],[size(leadAxonOrder,1)+0.5 ,size(leadAxonOrder,1)+0.5],'color','k','LineWidth',2);
line([size(leadAxonOrder,1)+0.5 ,size(leadAxonOrder,1)+0.5],[0,size(SaccadicCounts,2)],'color','k','LineWidth',2);
colorbar;
box on;
daspect([1,1,1]);




% for i = 1:size(SaccadicShortestPaths,2)
%     ii = SaccadicShortestPaths(i).ij(1);
%     jj = SaccadicShortestPaths(i).ij(2);
%    distMat(ii,jj) = SaccadicShortestPaths(i).dist;
% end
% distMat = distMat.*~eye(size(distMat));



%% Integrator Axons

UniqueDBXIntegratorAxons = unique(vertcat(DBX.Integrator));
DBXIntegratorMotorDiff = MotorDiff(UniqueDBXIntegratorAxons,df);
histogram(DBXIntegratorMotorDiff,10);

LeadLikeIntegratorAxons = UniqueDBXIntegratorAxons(DBXIntegratorMotorDiff>0);
LagLikeIntegratorAxons = UniqueDBXIntegratorAxons(DBXIntegratorMotorDiff<0);

[DBXIntegratorAxons,IntegratorCounts,IntegratorGraph,IntegratorShortestPaths] = PartnerConnectivity(UniqueDBXIntegratorAxons,UniqueDBXIntegratorAxons,'Integrator',df);
IntegratorCounts = IntegratorCounts.*~eye(size(IntegratorCounts));


leadAxonOrder = find(ismember(UniqueDBXIntegratorAxons,LeadLikeIntegratorAxons));
lagAxonOrder = find(ismember(UniqueDBXIntegratorAxons,LagLikeIntegratorAxons));

figure;
subplot(2,3,[1,4])
transform_swc_AV(LeadLikeIntegratorAxons,lightRed,[],true,false);
subplot(2,3,[2,5])
transform_swc_AV(LagLikeIntegratorAxons,lightBlue,[],true,false);
subplot(2,3,3)
imagesc(IntegratorCounts([leadAxonOrder;lagAxonOrder],[leadAxonOrder;lagAxonOrder]));
minMax = [min(IntegratorCounts(:)), max(IntegratorCounts(:))];
colorcet('L1','N',minMax(2)-minMax(1),'reverse',1);
line([0,size(IntegratorCounts,2)],[size(leadAxonOrder,1)+0.5 ,size(leadAxonOrder,1)+0.5],'color','k','LineWidth',2);
line([size(leadAxonOrder,1)+0.5 ,size(leadAxonOrder,1)+0.5],[0,size(IntegratorCounts,2)],'color','k','LineWidth',2);
colorbar;
box on;
daspect([1,1,1]);

% IntegratorLoopSize = zeros(size(IntegratorCounts));
% 
% for i = 1:size(IntegratorCounts,1)
%     ii = IntegratorShortestPaths(i).ij(1);
%     jj = IntegratorShortestPaths(i).ij(2);
%    IntegratorLoopSize(ii,jj) = IntegratorShortestPaths(i).size;
% end
% %distMat = distMat.*~eye(size(distMat));
% 
% subplot(2,3,6)
% histogram(IntegratorShortestPaths,'FaceColor','none','EdgeColor',colors(1,:),'LineWidth',4,'DisplayStyle','stairs');
% hold on;
% histogram(IntegratorShortestPaths,'FaceColor','none','EdgeColor',colors(2,:),'LineWidth',4,'DisplayStyle','stairs')
% box off;
% xlabel('shortest path nodes');
% ylabel('count');

%% Vestibular Axons

UniqueDBXVestibularAxons = unique(vertcat(DBX.Vestibular));
DBXVestibularMotorDiff = MotorDiff(UniqueDBXVestibularAxons,df);
histogram(DBXVestibularMotorDiff,10);

LeadLikeVestibularAxons = UniqueDBXVestibularAxons(DBXVestibularMotorDiff>0);
LagLikeVestibularAxons = UniqueDBXVestibularAxons(DBXVestibularMotorDiff<0);

[DBXVestibularAxons,VestibularCounts,VestibularGraph,VestibularShortestPaths] = PartnerConnectivity(UniqueDBXVestibularAxons,UniqueDBXVestibularAxons,'Vestibular',df);
VestibularCounts = VestibularCounts.*~eye(size(VestibularCounts));


leadAxonOrder = find(ismember(UniqueDBXVestibularAxons,LeadLikeVestibularAxons));
lagAxonOrder = find(ismember(UniqueDBXVestibularAxons,LagLikeVestibularAxons));

figure;
subplot(2,3,[1,4])
transform_swc_AV(LeadLikeVestibularAxons,lightRed,[],true,false);
subplot(2,3,[2,5])
transform_swc_AV(LagLikeVestibularAxons,lightBlue,[],true,false);
subplot(2,3,3)
imagesc(VestibularCounts([leadAxonOrder;lagAxonOrder],[leadAxonOrder;lagAxonOrder]));
minMax = [min(VestibularCounts(:)), max(VestibularCounts(:))];
colorcet('L1','N',minMax(2)-minMax(1),'reverse',1);
line([0,size(VestibularCounts,2)],[size(leadAxonOrder,1)+0.5 ,size(leadAxonOrder,1)+0.5],'color','k','LineWidth',2);
line([size(leadAxonOrder,1)+0.5 ,size(leadAxonOrder,1)+0.5],[0,size(VestibularCounts,2)],'color','k','LineWidth',2);
colorbar;
box on;
daspect([1,1,1]);

% IntegratorLoopSize = zeros(size(IntegratorCounts));
% 
% for i = 1:size(IntegratorCounts,1)
%     ii = IntegratorShortestPaths(i).ij(1);
%     jj = IntegratorShortestPaths(i).ij(2);
%    IntegratorLoopSize(ii,jj) = IntegratorShortestPaths(i).size;
% end
% %distMat = distMat.*~eye(size(distMat));
% 
% subplot(2,3,6)
% histogram(IntegratorShortestPaths,'FaceColor','none','EdgeColor',colors(1,:),'LineWidth',4,'DisplayStyle','stairs');
% hold on;
% histogram(IntegratorShortestPaths,'FaceColor','none','EdgeColor',colors(2,:),'LineWidth',4,'DisplayStyle','stairs')
% box off;
% xlabel('shortest path nodes');
% ylabel('count');