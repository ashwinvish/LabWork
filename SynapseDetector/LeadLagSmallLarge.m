clc;
clear;

if ismac
    addpath(genpath('/Users/ashwin/Documents/'));
    df = readtable('/Users/ashwin/Documents/SynapseDetector/11252018.csv');
else
    addpath(genpath('/usr/people/ashwinv/seungmount/research/Ashwin/Scripts'));
    df = readtable('/usr/people/ashwinv/seungmount/research/Ashwin/SynapseDetector/11252018.csv');
end

colorPallete = cbrewer('div','BrBG',5);

smallColor = colorPallete(1,:);
largeColor = colorPallete(5,:);

% load integrator neurons
load('LeadLikeDBX.mat');
load('LagLikeDBX.mat');
load('LeadLikeALX.mat');
load('LagLikeALX.mat');
load('leadNeurons.mat');
load('lagNeurons.mat');

% load Motor neurons
load('smallABDneurons.mat');
load('largeABDneurons.mat');
load('ABDr.mat');
load('ABDc.mat');
load('ABDIr.mat');
load('ABDIc.mat');


% load Saccadic neurons

load('LeadFunctionalSaccadicAxons.mat');
load('LagFunctionalSaccadicAxons.mat');
load('LeadLikeALXSaccadicAxons.mat');
load('LeadLikeDBXSaccadicAxons.mat');


colorPallete = colorcet('R3','N',4);

% load saccadic cell IDs
triangle = [76749 76752 77374 77446 77455 78679 ];

sparseSaccadic = [76540 76622 76626 76697 76748 76750 76751 77122 77151 ...
    77238 77239 77240 77241 77437 77645 77708 77740 77826 78351 78545 78558 ...
    78572 78641 76611 77630 77374 80763 80850 80821 80801 80743 81007 80974 ...
    81002 80216 ];

bushySaccadicMedial = [76618 76625 76627 77132 77162 77163 77329 77434 77447 77460 ...
    77467 77797 77805 77848 78357 78358 79054 79059 78544 78650 77390 77621 ...
    77636 77651 77656 80995 ];

bushySaccadicLateral = [79067 79058 79069 79072 79055 79074 79062 79077 79085 ...
    79080 79086 79064 78576 78540 77146 78646 78542 78546 78541 78543 80728]; 

putativeBushySaccadic = [77342 77352 77336 77354 77373 78346 77433 77435 77453 77461];

lateralDSaccadic = [76629 76667 80185 76691 76692 77357 77389 77689 77806 ...
                    77816 78601 77667 77684 80542 ];

lateralVSaccadic = [80163 80167 80177 80179 76688 80204 80206 80210 76682];

unknownSaccadic = [78583 78649 77670 80217 80746 80679 80804 80757 80647 ...
    80947 80939 81027 80943 80315 ];

allSaccadic = [triangle,sparseSaccadic,bushySaccadicMedial,bushySaccadicLateral,...
    putativeBushySaccadic,lateralDSaccadic,lateralVSaccadic,unknownSaccadic];

functionalBurstNeurons = [leadNeurons';LeadLikeDBX;LeadLikeALX];
functionalTonicNeurons = [lagNeurons';LagLikeDBX;LagLikeALX];

% Motor Neuron Cell IDs

ABDr_CellIDs = [77648, 77710, 77300, 77705, 77305, 77301, 77709, 77672, 77302];
ABDc_CellIDs = [77154, 77646, 77682 ,77628 ,77295 , 77652 ,77292 ,77688 ,77654 ,77658 ,77657 ,77662, 77296];
ABDIr_CellIDs = [77631, 77150, 77618, 77886, 78547, 77158, 78556, 78552, 77665, 77668, 77634, 78553];
ABDIc_CellIDs = [77148, 77625, 77641, 77692, 77144, 77643, 77640, 79051, 79066, 78574];

AllABD = [ABDr_CellIDs,ABDc_CellIDs,ABDIr_CellIDs,ABDIc_CellIDs];
Allmotor = [ABDr_CellIDs,ABDc_CellIDs];

%% classify all saccadic axons onto small and Lagre ABd neurons

for i = 1:numel(allSaccadic)
    Saccade(i).cellID = allSaccadic(i);
    [A,~] = SynapticPartners(allSaccadic(i),2,df);
    A = A(A<1e5);
    Saccade(i).isPartnerLead = ismember(A,functionalBurstNeurons);
    Saccade(i).isPartnerLag = ismember(A,functionalTonicNeurons);
    Saccade(i).Lead = A(Saccade(i).isPartnerLead);
    Saccade(i).Lag = A(Saccade(i).isPartnerLag);
    Saccade(i).MotorDist = isMotor(allSaccadic(i),df);
    Saccade(i).isPartnersmall = ismember(A,smallABDCellIDs);
    Saccade(i).isPartnerLarge = ismember(A,largeABDCellIDs);
    Saccade(i).SmallABD = A(Saccade(i).isPartnersmall);
    Saccade(i).LargeABD = A(Saccade(i).isPartnerLarge);
    clear A;
end


%% plots

for i = 1:numel(allSaccadic)
    axon.ID(i) = allSaccadic(i);
    axon.Conn(1,1,i) = numel(Saccade(i).SmallABD) + numel(Saccade(i).Lead);
    axon.Conn(2,1,i) = numel(Saccade(i).LargeABD) + numel(Saccade(i).Lead);
    axon.Conn(1,2,i) = numel(Saccade(i).SmallABD) + numel(Saccade(i).Lag);
    axon.Conn(2,2,i) = numel(Saccade(i).LargeABD) + numel(Saccade(i).Lag);
    axon.mean(i) = mean(mean(axon.Conn(:,:,i)));
    if ~isempty(Saccade(i).SmallABD)
        axon.SmallLarge(i) = numel(Saccade(i).SmallABD) / (numel(Saccade(i).SmallABD)+numel(Saccade(i).LargeABD));
    else
        axon.SmallLarge(i) = NaN;
    end
    
     if ~isempty(Saccade(i).LargeABD)
        axon.LargeSmall(i) = numel(Saccade(i).LargeABD) / (numel(Saccade(i).SmallABD)+numel(Saccade(i).LargeABD));
    else
        axon.LargeSmall(i) = NaN;
     end
    
    if ~isempty(Saccade(i).Lead)
        axon.BurstTonic(i) = numel(Saccade(i).Lead) / (numel(Saccade(i).Lead)+ numel(Saccade(i).Lag));
    else
        axon.BurstTonic(i) = NaN;
    end
end

figure;
subplot(4,4,1)
plot(repmat([1;2],1,135),[axon.BurstTonic;axon.SmallLarge],'-ko');
set(gca,'XLim',[0,3],'XTick',[1,2],'XTickLabel',[{'B/T'},{'S/L'}]);
ylabel('ns/(ns+nl) ; nb/(nb+nt)')
line([0,3],[0.5,0.5],'color','r','LineStyle','--');
axis square;
box off;

subplot(4,4,2)
plot(repmat([1;2],1,135),[axon.BurstTonic;axon.LargeSmall],'-ko');
set(gca,'XLim',[0,3],'XTick',[1,2],'XTickLabel',[{'B/T'},{'L/S'}]);
line([0,3],[0.5,0.5],'color','r','LineStyle','--');
axis square;
box off;

% 

for i = 1:numel(allSaccadic)
    if sum(Saccade(i).Lead) > sum(Saccade(i).Lag)
        LeadSaccade(i) = allSaccadic(i);
    else
        sum(Saccade(i).Lead) < sum(Saccade(i).Lag)
        LagSaccade(i) = allSaccadic(i);
    end
end
LeadSaccade = LeadSaccade(LeadSaccade~=0);
LagSaccade = LagSaccade(LagSaccade~=0);

save('LeadLikeSaccadicAxons.mat','LeadSaccade');
save('LagLikeSaccadicAxons.mat','LagSaccade');

figure;
% transform_swc_AV(LeadSaccade,'r',[],true,false);
% transform_swc_AV(LagSaccade,'b',[],false,false);

%% classify all Saccadic as Lead Lag

for i = 1:numel(smallABDCellIDs)
    if ismember(smallABDCellIDs(i), ABDr_CellIDs)
        index = find(smallABDCellIDs(i) == ABDr_CellIDs);
        motorNeuron = ABDr(index);
    else
        index = find(smallABDCellIDs(i) == ABDc_CellIDs);
        motorNeuron = ABDc(index);
    end
    smallSaccadicAxons(i).smallABDID = smallABDCellIDs(i);
    smallSaccadicAxons(i).SaccadicIDs = motorNeuron.Saccadic;
    smallSaccadicAxons(i).PathLength = motorNeuron.PathLength;
    
    smallSaccadicAxons(i).LeadPathLength = [];
    smallSaccadicAxons(i).LagPathLength = [];
    smallSaccadicAxons(i).LeadPathLengthNorm = [];
    smallSaccadicAxons(i).LagPathLengthNorm = [];
    
    for j = 1:numel(smallSaccadicAxons(i).SaccadicIDs)
        if ismember(smallSaccadicAxons(i).SaccadicIDs(j),LeadSaccade)
            ixLead = find(smallSaccadicAxons(i).SaccadicIDs(j) == motorNeuron.Inputs);
            smallSaccadicAxons(i).LeadPathLength = [smallSaccadicAxons(i).LeadPathLength; motorNeuron.PathLength(ixLead)];
            smallSaccadicAxons(i).LeadPathLengthNorm = histcounts(motorNeuron.PathLength(ixLead)./max(Pvec_tree(motorNeuron.Tree{1})),0:0.1:1,'Normalization','probability');
            %smallSaccadicAxons(i).UniqueAxons = 
            clear ixLead;
        elseif ismember(smallSaccadicAxons(i).SaccadicIDs(j),LagSaccade)
            ixLag = find(smallSaccadicAxons(i).SaccadicIDs(j) == motorNeuron.Inputs);
            smallSaccadicAxons(i).LagPathLength = [smallSaccadicAxons(i).LagPathLength; motorNeuron.PathLength(ixLag)];
            smallSaccadicAxons(i).LagPathLengthNorm = histcounts(motorNeuron.PathLength(ixLag)./max(Pvec_tree(motorNeuron.Tree{1})),0:0.1:1,'Normalization','probability');
            clear ixLag;
        end
    end
    smallSaccadicAxons(i).LeadLagRatio = size(smallSaccadicAxons(i).LeadPathLength,1) / size(smallSaccadicAxons(i).LagPathLength,1);
    clear motorNeuron;
end


for i = 1:numel(largeABDCellIDs)
    if ismember(largeABDCellIDs(i), ABDr_CellIDs)
        index = find(largeABDCellIDs(i) == ABDr_CellIDs);
        motorNeuron = ABDr(index);
    else
        index = find(largeABDCellIDs(i) == ABDc_CellIDs);
        motorNeuron = ABDc(index);
    end
    largeSaccadicAxons(i).smallABDID = largeABDCellIDs(i);
    largeSaccadicAxons(i).SaccadicIDs = motorNeuron.Saccadic;
    largeSaccadicAxons(i).PathLength = motorNeuron.PathLength;
    
    largeSaccadicAxons(i).LeadPathLength = [];
    largeSaccadicAxons(i).LagPathLength = [];
    largeSaccadicAxons(i).LeadPathLengthNorm = [];
    largeSaccadicAxons(i).LagPathLengthNorm = [];
    
    for j = 1:numel(largeSaccadicAxons(i).SaccadicIDs)
        if ismember(largeSaccadicAxons(i).SaccadicIDs(j),LeadSaccade)
            ixLead = find(largeSaccadicAxons(i).SaccadicIDs(j) == motorNeuron.Inputs);
            largeSaccadicAxons(i).LeadPathLength = [largeSaccadicAxons(i).LeadPathLength; motorNeuron.PathLength(ixLead)];
            largeSaccadicAxons(i).LeadPathLengthNorm = histcounts(motorNeuron.PathLength(ixLead)./max(Pvec_tree(motorNeuron.Tree{1})),0:0.1:1,'Normalization','probability');
            clear ixLead;
        elseif ismember(largeSaccadicAxons(i).SaccadicIDs(j),LagSaccade)
            ixLag = find(largeSaccadicAxons(i).SaccadicIDs(j) == motorNeuron.Inputs);
            largeSaccadicAxons(i).LagPathLength = [largeSaccadicAxons(i).LagPathLength; motorNeuron.PathLength(ixLag)];
            largeSaccadicAxons(i).LagPathLengthNorm = histcounts(motorNeuron.PathLength(ixLag)./max(Pvec_tree(motorNeuron.Tree{1})),0:0.1:1,'Normalization','probability');
            clear ixLag;
        end
    end
    largeSaccadicAxons(i).LeadLagRatio = size(largeSaccadicAxons(i).LeadPathLength ,1) / size(largeSaccadicAxons(i).LagPathLength ,1);
    clear motorNeuron;
end



%% plots
figure;
% smallABD
subplot(4,4,1)
histogram(vertcat(smallSaccadicAxons.LeadPathLength),20,'Normalization','probability','EdgeColor',smallColor,'DisplayStyle','stairs','LineWidth',2);
hold on;
histogram(vertcat(largeSaccadicAxons.LeadPathLength),20,'Normalization','probability','EdgeColor',largeColor,'DisplayStyle','stairs','LineWidth',2);
box off;
axis square;
title('Burst like');
legend({'Small (mif)','Large (SIF)'},'Location','bestoutside');

subplot(4,4,2)
histogram(vertcat(smallSaccadicAxons.LagPathLength),20,'Normalization','probability','EdgeColor',smallColor,'DisplayStyle','stairs','LineWidth',2);
hold on;
histogram(vertcat(largeSaccadicAxons.LagPathLength),20,'Normalization','probability','EdgeColor',largeColor,'DisplayStyle','stairs','LineWidth',2);
box off;
title('Tonic like');
axis square;

% subplot(4,4,3)
% histogram(vertcat(smallSaccadicAxons.LeadPathLengthNorm),20,'Normalization','probability','EdgeColor',smallColor,'DisplayStyle','stairs','LineWidth',2);
% hold on;
% histogram(vertcat(largeSaccadicAxons.LeadPathLengthNorm),20,'Normalization','probability','EdgeColor',largeColor,'DisplayStyle','stairs','LineWidth',2);
% box off;
% axis square;
% title('Burst like');
% legend({'Small (mif)','Large (SIF)'},'Location','bestoutside');
% 
% subplot(4,4,4)
% histogram(vertcat(smallSaccadicAxons.LagPathLengthNorm),20,'Normalization','probability','EdgeColor',smallColor,'DisplayStyle','stairs','LineWidth',2);
% hold on;
% histogram(vertcat(largeSaccadicAxons.LagPathLengthNorm),20,'Normalization','probability','EdgeColor',largeColor,'DisplayStyle','stairs','LineWidth',2);
% box off;
% title('Tonic like');
% axis square;

% Errorbar plots


subplot(4,4,5)
h1 = shadedErrorBar([0.1:0.1:1],mean(vertcat(smallSaccadicAxons.LeadPathLengthNorm)),std(vertcat(smallSaccadicAxons.LeadPathLengthNorm))./sqrt(numel(smallABDCellIDs)),...
    'lineProps',{'Color',smallColor,'LineWidth',2});
hold on;
h2 = shadedErrorBar([0.1:0.1:1],mean(vertcat(largeSaccadicAxons.LeadPathLengthNorm)),std(vertcat(largeSaccadicAxons.LeadPathLengthNorm))./sqrt(numel(largeABDCellIDs)),...
    'lineProps',{'Color',largeColor,'LineWidth',2});box off;
axis square;
title('Burst like');
legend([h1.patch.Parent], {'Small (MIF)','Large (SIF)'},'Location','bestoutside');

subplot(4,4,6)
h1 = shadedErrorBar([0.1:0.1:1],mean(vertcat(smallSaccadicAxons.LagPathLengthNorm)),std(vertcat(smallSaccadicAxons.LagPathLengthNorm))./sqrt(numel(smallABDCellIDs)),...
    'lineProps',{'Color',smallColor,'LineWidth',2});
hold on;
h2 = shadedErrorBar([0.1:0.1:1],mean(vertcat(largeSaccadicAxons.LagPathLengthNorm)),std(vertcat(largeSaccadicAxons.LagPathLengthNorm))./sqrt(numel(largeABDCellIDs)),...
    'lineProps',{'Color',largeColor,'LineWidth',2});box off;
axis square;
title('tonic like');

