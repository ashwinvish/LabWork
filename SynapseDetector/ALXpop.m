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

confirmedALX = [76181 76201 76187 76184 76192 76197 ];
putativeALX = [76657 76623 76701 76661 76673 76690 76677 76624 77327 77580 ...
    77341 77344 77868 77872 77349 76634 77592 77578 77607 78629 77581 77591 ...
    77602 77821 77844 77845 77822 78406 78667 79022 79033 78404 78441 78421 ...
    79044 79046 79048 80221 78853 79017 79852 78451 79042 80596 80606 78911 ...
    79746 80271 79720 79976 77586 77369 78633 80750 77142 79060 78453 80885];

allALX = [confirmedALX,putativeALX];

for i = 1:numel(allALX)
    ALX(i) = InputsByClass(allALX(i),df);
end

uniqeALXSaccadicAxons = unique(vertcat(ALX.Saccadic));
ALXSaccadicMotorDist = MotorDiff(uniqeALXSaccadicAxons,df);

%histogram(ALXSaccadicMotorDist,10);

LeadLikeALXSaccadicAxons = uniqeALXSaccadicAxons(ALXSaccadicMotorDist>0);
LagLikeALXSaccadicAxons = uniqeALXSaccadicAxons(ALXSaccadicMotorDist<0);
subplot(2,3,[1,4])
transform_swc_AV(LeadLikeALXSaccadicAxons,colors(1,:),[],true,false);
subplot(2,3,[2,5])
transform_swc_AV(LagLikeALXSaccadicAxons,colors(2,:),[],true,false);

%% partners of lead and Lag axons and how recurrently they.

[LeadLikeAxons,countsLead,Leadgraph,LeadShortestPaths] = PartnerConnectivity(LeadLikeALXSaccadicAxons,uniqeALXSaccadicAxons,df);
[LagLikeAxons,countsLag,Laggraph,LagShortestPaths] = PartnerConnectivity(LagLikeALXSaccadicAxons,uniqeALXSaccadicAxons,df);

countsLeadLag = vertcat(countsLead,countsLag);

subplot(2,3,3)
imagesc(countsLeadLag);
minMax = [min(countsLeadLag(:)), max(countsLeadLag(:))];
colorcet('L1','N',minMax(2)-minMax(1),'reverse',1);
line([0,size(countsLeadLag,2)],[size(LeadLikeAxons,2)+0.5 ,size(LeadLikeAxons,2)+0.5],'color','k','LineWidth',2);
line([size(LeadLikeAxons,2)+0.5 ,size(LeadLikeAxons,2)+0.5],[0,size(countsLeadLag,2)],'color','k','LineWidth',2);
colorbar;
box on;
daspect([1,1,1]);

subplot(2,3,6)
histogram(LeadShortestPaths,'FaceColor','none','EdgeColor',colors(1,:),'LineWidth',4,'DisplayStyle','stairs');
hold on;
histogram(LagShortestPaths,'FaceColor','none','EdgeColor',colors(2,:),'LineWidth',4,'DisplayStyle','stairs')
box off;
xlabel('shortest path nodes');
ylabel('count');

%%






