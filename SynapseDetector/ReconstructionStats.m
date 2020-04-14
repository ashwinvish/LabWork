
% reconstructed fraction

addpath(genpath('/Users/ashwin/Documents/'));
if ismac
    addpath(genpath('/Users/ashwin/Documents/'));
    df = readtable('/Users/ashwin/Google Drive/Zfish/SynapseDetector/04152019.csv');
else
    addpath(genpath('/usr/people/ashwinv/seungmount/research/Ashwin/Scripts'));
    df = readtable('/usr/people/ashwinv/seungmount/research/Ashwin/SynapseDetector/11252018.csv');
end



load SaccadicProjectingToABDexclusively.mat
load SaccadicProjectingToABDiexclusively.mat

r456ABD = reconstructedFraction(SaccadicProjectingToABDexclusively.cellID,df);
r456ABDi = reconstructedFraction(SaccadicProjectingToABDiexclusively.cellID,df);

load ipsir78Int.mat
load contrar78Int.mat

r78ipsi = reconstructedFraction(ipsir78Int,df);
r78contra = reconstructedFraction(contrar78Int,df);

load ABDr.mat
load ABDc.mat

ABDrRecon = reconstructedFraction(vertcat(ABDr.cellID),df);
ABDcRecon = reconstructedFraction(vertcat(ABDc.cellID),df);


figure;
subplot(4,4,1);
plot(ones(numel(r456ABD),1),r456ABD,'o');
hold on
plot(2*ones(numel(r456ABDi),1),r456ABDi,'o');
plot(3*ones(numel(r78ipsi),1),r78ipsi,'o');
plot(4*ones(numel(r78contra),1),r78contra,'o');
plot(5*ones(numel(ABDrRecon),1),ABDrRecon,'o');
plot(6*ones(numel(ABDcRecon),1),ABDcRecon,'o');
line([0,6],[0.8,0.8])
legend({'r456_M','r456_I','r78_ipsi','r78_contra','ABDr','ADBc'})
box off;
axis square;


%%

confirmedALX = [76181 76201 76187 76184 76192 76197 ];
putativeALX = [76657 76623 76701 76661 76673 76690 76677 76624 77327 77580 ...
    77341 77344 77868 77872 77349 76634 77592 77578 77607 78629 77581 77591 ...
    77602 77821 77844 77845 77822 78406 78667 79022 79033 78404 78441 78421 ...
    79044 79046 79048 80221 78853 79017 79852 78451 79042 80596 80606 78911 ...
    79746 80271 79720 79976 77586 77369 78633 80750 77142 79060 78453 80885 ...
    81423 81661 81683 81792];

allPutativeALX = putativeALX;
putativeALX = putativeALX(logical(isPostSynapseIntegrator(putativeALX,df))|logical(isPreSynapseIntegrator(putativeALX,df)));
motorOut = isMotor(putativeALX',df);
allPutativeALX_beforeMotorFilter = putativeALX;
putativeALX = putativeALX(sum(motorOut(:,2:end),2)>0);
ALX.cellIDs = [confirmedALX,putativeALX];
ALX.cellIDs = ALX.cellIDs(isExistReRoot(ALX.cellIDs));

ALXexcluded = setdiff(allPutativeALX_beforeMotorFilter,putativeALX)';

% transform_swc_AV(ALX.cellIDs,ALXcolor,[],true,true);
% transform_swc_AV(ALXexcluded,[254,246,91]./255,[],false,false);


%%
load SaccadicProjectingToABDexclusively.mat
load SaccadicProjectingToABDiexclusively.mat
load SaccadicProjectingToABD_ABDi.mat

bushySaccadicLateral = [79067 79058 79069 79072 79055 79074 79062 79077 79085 ...
    79080 79086 79064 78576 78540 77146 78646 78542 78546 78541 78543 80728]; 


% transform_swc_AV(SaccadicProjectingToABDiexclusively.cellID,SaccABDicolor,[],true,true);
% transform_swc_AV(SaccadicProjectingToABDexclusively.cellID,SaccABDcolor,[],false,false);
% transform_swc_AV(bushySaccadicLateral,[254,246,91]./255,[],false,false);

%%

load AllCells.mat

% for i = 1:numel(AllCells)
%     if isExistReRoot(AllCells(i))
%         if ~isContra(AllCells(i))
%         allCells(i) = isRhombomere(AllCells(i));
%         end
%     end
% end
% IpsiNeurons = vertcat(allCells.cellID);
% numberOfIpsineurons = size(IpsiNeurons,1);
% preInteg = isPreSynapseIntegrator(IpsiNeurons,df);
% postInteg = isPostSynapseIntegrator(IpsiNeurons,df);
% 
% recurrentNeurons = IpsiNeurons(find(preInteg+postInteg == 2));

allApprovedIntegrators = [ALX.cellIDs'; ...
    SaccadicProjectingToABDexclusively.cellID';...
    SaccadicProjectingToABDiexclusively.cellID';...
    SaccadicProjectingToABD_ABDi.cellID'];


setdiff(recurrentNeurons,allApprovedIntegrators)

temp = [76613,77156,76554,76871,83327,76555,76902,81004,76872,77157,82807,79373,79948,...
    78083,79962,76555,79948,76887,...
    77889,77157,77156,78566,79211,76614];





