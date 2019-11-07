 % input distribution
 
 % Classes by anatomy
clear;
addpath(genpath('/Users/ashwin/Documents/'));
colorSchemes;

colors = cbrewer('qual','Paired',10);

temp1 = colorcet('CBTL1','N',5); %  linear-tritanopic_krjcw_5-98_c46_n256
temp2 = colorcet('CBL2','N',5);  %  linear-protanopic-deuteranopic_kbw_5-98_c40_n256

leadColor = temp1(3,:); % dark red, CB safe
lagColor = temp2(3,:); % dark blue, CB safe

PartnerColors = colorcet('CBTD1','N',10);
lightRed = PartnerColors(10,:);
lightBlue = PartnerColors(1,:);

colorPallete = colorcet('D2','N',5);
colorSchemes
lightGreen = colorPallete(1,:);
lightMagenta = colorPallete(5,:);


startup

if ismac
    addpath(genpath('/Users/ashwin/Documents/'));
    df = readtable('/Users/ashwin/Google Drive/Zfish/SynapseDetector/04152019.csv');
else
    addpath(genpath('/usr/people/ashwinv/seungmount/research/Ashwin/Scripts'));
    df = readtable('/usr/people/ashwinv/seungmount/research/Ashwin/SynapseDetector/11252018.csv');
end

%% Load Saccadic Neurons

triangle = [76749 76752 77446 77455 78679];

sparseSaccadic = [76540,76622,76626,76697,76748,76750,76751,77122,77151,77238,...
    77239,77240,77241,77437,77645,77708,77740,77826,78351,78558,78572,78641,76611,...
    77630,77374,80763,80850,80821,80801,80743,81007,80974,81002,80216,80681,81161,81145,77304,81793];

bushySaccadicMedial = [76618,76625,76627,77132,77162,77163,77329,77434,77447,...
    77460,77467,77797,77805,77848,78357,78358,79054,79059,78544,78650,77390,77621,...
    77636,77651,77656,80995,81406,80629,81559,81363,81293,81338,80956,81395,81407,...
    81503,81312,81297,81417,79078,81295,81317,79085,81336,80285,];

%bushySaccadicLateral = [79067 79058 79069 79072 79055 79074 79062 79077 79085 ...
%    79080 79086 79064 78576 78540 77146 78646 78542 78546 78541 78543 80728]; 

putativeBushySaccadic = [77342 77352 77336 77354 77373 78346 77433 77435 77453 77461];

lateralDSaccadic = [76629 76667 80185 76691 76692 77357 77389 77689 77806 ...
                    77816 78601 77667 77684 80542 ];

lateralVSaccadic = [80163 80167 80177 80179 76688 80204 80206 80210 76682];

unknownSaccadic = [78583 78649 77670 80217 80746 80679 80804 80757 80647 ...
    80947 80939 81027 80943 80315 ];

%SaccadicAxons = [triangle,sparseSaccadic,bushySaccadicMedial,bushySaccadicLateral,...
%    putativeBushySaccadic,lateralDSaccadic,lateralVSaccadic,unknownSaccadic];

SaccadicAxons = [triangle,sparseSaccadic,bushySaccadicMedial,...
    putativeBushySaccadic,lateralDSaccadic,unknownSaccadic];

%IBNordered = [77125 77128 77231 77941 77940 77157 77247 77153 79053 77135 77941 78550];
%IBNrem = [77137 77942 78557 78567 78685 79083 79084 80971];
IBNall = [ 77125 77128 77135 77153 77231 77247 78685 77941 77137 79053 77940 77942 80971 77157 78550 79084 78557 79083 78567];

confirmedALX = [76181 76201 76187 76184 76192 76197 ];
putativeALX = [76657 76623 76701 76661 76673 76690 76677 76624 77327 77580 ...
    77341 77344 77868 77872 77349 76634 77592 77578 77607 78629 77581 77591 ...
    77602 77821 77844 77845 77822 78406 78667 79022 79033 78404 78441 78421 ...
    79044 79046 79048 80221 78853 79017 79852 78451 79042 80596 80606 78911 ...
    79746 80271 79720 79976 77586 77369 78633 80750 77142 79060 78453 80885 ...
    81423 81661 81683 81792 76389,76873,77935,79743,81374];
putativeALX = putativeALX(logical(isPostSynapseIntegrator(putativeALX,df))|logical(isPreSynapseIntegrator(putativeALX,df)));
motorOut = isMotor(putativeALX',df);
putativeALX = putativeALX(sum(motorOut(:,2:end),2)>0);

ALX.cellIDs = [confirmedALX,putativeALX];
ALX.cellIDs = ALX.cellIDs(isExistReRoot(ALX.cellIDs));

confirmedDBX = [76199 76200 76182 76183 76185 76186 76189 76191 76188 ];
putativeDBX = [77582 77588 77595 77597 77598 77599 77605 79040 76399 76828 ...
    76829 78435 76826 76874 76289 76542 76832 76836 76838 76877 76883 76884 ...
    78400 78429 76383 76867 77594 76830 76876 78440 78438 78447 79045 78434 ...
    78450 76878 77683 77450 78448 ];

putativeDBX = putativeDBX(logical(isPreSynapseIntegrator(putativeDBX,df)));


confirmedBARHL = [76198 76190 76193 76194 76195 76196 ];
putativeBARHL = [78452 80224 78391];

allALX = [confirmedALX,putativeALX];
allDBX = [confirmedDBX,putativeDBX];

%%
ABDr_CellIDs = [82140,82145,82143,77648,82146,77302,77301, 82194,77705,77710,77305,77672,77300,77709];
ABDc_CellIDs = [82213,77654,77646,82212,77295,81172,77657,82197,77682,77154, 77652,77658,77628,77292,77688,82195,77296];
ABDIr_CellIDs = [78553,77150,77886,78547,77631,77158,77665,78552,77668,77618,77634];
ABDIc_CellIDs = [78574,79051,77148,77625,77144,77692,77641,79066,77640,77643];


load ABDVols.mat
load ABDr.mat
load ABDc.mat
load ABDIr.mat
load ABDIc.mat

%%
 
 for i = 1:numel(SaccadicAxons)
     if ~isExistReRoot(SaccadicAxons(i))==0
     Sacc(i) = InputsByClass(SaccadicAxons(i),df,1);
     end  
 end
 
 
 load SaccadicProjectingToABDexclusively.mat
 load SaccadicProjectingToABDiexclusively.mat
 load ABDPutativeSaccadic.mat
 load ABDiPutativeSaccadic.mat
 
 %  
%% Cumulative Plots of inputs

% find cellIDs in Sacc;

SaccadicProjectingToABDexclusively.IDlocs = find(ismember(SaccadicAxons,SaccadicProjectingToABDexclusively.cellID));


[SaccadicProjectingToABDexclusively.VestibularPathLength, SaccadicProjectingToABDexclusively.VestibularGradient] = ...
    getABDgradient(Sacc(SaccadicProjectingToABDexclusively.IDlocs),unique(vertcat(Sacc.Vestibular)),true,false);

[SaccadicProjectingToABDexclusively.r78PathLength, SaccadicProjectingToABDexclusively.r78Gradient] = ...
    getABDgradient(Sacc(SaccadicProjectingToABDexclusively.IDlocs),ALX.cellIDs',true,false);

[SaccadicProjectingToABDexclusively.r56PathLength, SaccadicProjectingToABDexclusively.r56Gradient] = ...
    getABDgradient(Sacc(SaccadicProjectingToABDexclusively.IDlocs),lateralVSaccadic,true,false);

[SaccadicProjectingToABDexclusively.r456PathLength, SaccadicProjectingToABDexclusively.r456Gradient] = ...
    getABDgradient(Sacc(SaccadicProjectingToABDexclusively.IDlocs),SaccadicProjectingToABDexclusively.cellID',true,false);

[SaccadicProjectingToABDexclusively.r23PathLength, SaccadicProjectingToABDexclusively.r23Gradient] = ...
    getABDgradient(Sacc(SaccadicProjectingToABDexclusively.IDlocs),[ABDPutativeSaccadic.cellIDs;ABDiPutativeSaccadic.cellIDs],true,false);


%ABDi

SaccadicProjectingToABDiexclusively.IDlocs = find(ismember(SaccadicAxons,SaccadicProjectingToABDiexclusively.cellID));


[SaccadicProjectingToABDiexclusively.VestibularPathLength, SaccadicProjectingToABDiexclusively.VestibularGradient] = ...
    getABDgradient(Sacc(SaccadicProjectingToABDiexclusively.IDlocs),unique(vertcat(Sacc.Vestibular)),true,false);

[SaccadicProjectingToABDiexclusively.r78PathLength, SaccadicProjectingToABDiexclusively.r78Gradient] = ...
    getABDgradient(Sacc(SaccadicProjectingToABDiexclusively.IDlocs),ALX.cellIDs',true,false);

[SaccadicProjectingToABDiexclusively.r56PathLength, SaccadicProjectingToABDiexclusively.r56Gradient] = ...
    getABDgradient(Sacc(SaccadicProjectingToABDiexclusively.IDlocs),lateralVSaccadic,true,false);

[SaccadicProjectingToABDiexclusively.r456PathLength, SaccadicProjectingToABDiexclusively.r456Gradient] = ...
    getABDgradient(Sacc(SaccadicProjectingToABDiexclusively.IDlocs),SaccadicProjectingToABDiexclusively.cellID',true,false);

[SaccadicProjectingToABDiexclusively.r23PathLength, SaccadicProjectingToABDiexclusively.r23Gradient] = ...
    getABDgradient(Sacc(SaccadicProjectingToABDiexclusively.IDlocs),[ABDPutativeSaccadic.cellIDs;ABDiPutativeSaccadic.cellIDs],true,false);




figure;
subplot(4,4,1)

histogram(cell2mat(SaccadicProjectingToABDexclusively.r23PathLength),20,'Normalization','cdf','DisplayStyle','stairs','EdgeColor','k','LineWidth',2);
hold on;
histogram(cell2mat(SaccadicProjectingToABDexclusively.VestibularPathLength),20,'Normalization','cdf','DisplayStyle','stairs','EdgeColor','r','LineWidth',2);
histogram([cell2mat(SaccadicProjectingToABDexclusively.r78PathLength)';cell2mat(SaccadicProjectingToABDexclusively.r56PathLength)';cell2mat(SaccadicProjectingToABDexclusively.r456PathLength)'],...
    20,'Normalization','cdf','DisplayStyle','stairs','EdgeColor','b','LineWidth',2);
xlabel('Norm. pathlength');
ylabel('Cumulative Count');
axis square;
box off;
offsetAxes(gca);


subplot(4,4,2)

histogram(cell2mat(SaccadicProjectingToABDiexclusively.r23PathLength),20,'Normalization','cdf','DisplayStyle','stairs','EdgeColor','k','LineWidth',2);
hold on;
histogram(cell2mat(SaccadicProjectingToABDiexclusively.VestibularPathLength),20,'Normalization','cdf','DisplayStyle','stairs','EdgeColor','r','LineWidth',2);
histogram([cell2mat(SaccadicProjectingToABDiexclusively.r78PathLength)';cell2mat(SaccadicProjectingToABDiexclusively.r56PathLength)';cell2mat(SaccadicProjectingToABDiexclusively.r456PathLength)'],...
    20,'Normalization','cdf','DisplayStyle','stairs','EdgeColor','b','LineWidth',2);
xlabel('Norm. pathlength');
ylabel('Cumulative Count');
axis square;
box off;
offsetAxes(gca);

 