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

% ABDvols = sortrows(ABDvols,2);
% ABDvols(:,2) = ABDvols(:,2)./1e9;
% 
% for i = 1:size(ABDvols,1)
%     ABDvols(i,3) = size(SynapticPartners(ABDvols(i,1),1,df),1);
%     if (isExistReRoot(ABDvols(i,1)) == 1)
%         tempTree = SwctoZbrian(ABDvols(i,1));
%         ABDvols(i,4) = max(Pvec_tree(tempTree{1}));
%         clear tempTree;
%         ABDvols(i,5:7) = getOrigin(ABDvols(i,1));
%     else
%         ABDvols(i,4)= NaN;  
%     end
%     %ABDvols(i,5:8) = getOrigin(ABDvols(i,1));
% end
% 
% ABDr_vols = ABDvols(find(ismember(ABDvols(:,1),ABDr_CellIDs)),2);
% ABDr_TotalSynapses = ABDvols(find(ismember(ABDvols(:,1),ABDr_CellIDs)),3);
% 
% ABDc_vols = ABDvols(find(ismember(ABDvols(:,1),ABDc_CellIDs)),2);
% ABDc_TotalSynapses = ABDvols(find(ismember(ABDvols(:,1),ABDc_CellIDs)),3);
% 
% ABDIr_vols = ABDvols(find(ismember(ABDvols(:,1),ABDIr_CellIDs)),2);
% ABDIr_TotalSynapses = ABDvols(find(ismember(ABDvols(:,1),ABDIr_CellIDs)),3);
% 
% ABDIc_vols = ABDvols(find(ismember(ABDvols(:,1),ABDIc_CellIDs)),2);
% ABDIc_TotalSynapses = ABDvols(find(ismember(ABDvols(:,1),ABDIc_CellIDs)),3);
% 
% figure;
% subplot(4,4,1)
% scatter(ABDvols(:,2),ABDvols(:,3),'ko');
% xlabel('ABD vol (um^3)');
% ylabel('Total synapses');
% axis square
% subplot(4,4,3)
% scatter(ABDvols(:,2),ABDvols(:,4),'ko');
% xlabel('ABD vol (um^3)');
% ylabel('ABD path length (um)');
% axis square
% subplot(4,4,6)
% scatter(ABDvols(:,3),ABDvols(:,4),'ko');
% scatter(ABDvols(:,4),ABDvols(:,3),'ko');
% xlabel('path length (um)');
% ylabel('Total synapses');
% axis square
% subplot(4,4,8)
% scatter3(ABDvols(:,2),ABDvols(:,4),ABDvols(:,3),'ko');
% xlabel('ABD vol (um^3)');
% ylabel('Pathlength (um)');
% zlabel('Total Synapses');
% axis square
% box on

%%
 for i = 1:numel(SaccadicAxons)
     if ~isExistReRoot(SaccadicAxons(i))==0
     Sacc(i) = InputsByClass(SaccadicAxons(i),df,2);
     %Sacc(i).Outputs = SynapticPartners(SaccadicAxons(i),2,df);
     end  
 end
 
% to get the input distribution onto the neurons


%%

Saccadic.MotorDistributionAllCounts = isMotor(SaccadicAxons',df);
%noSynapsesOnMotor = find(sum(SaccadicMotorDistribution.allCounts(:,2:5),2)==0);
%SaccadicMotorDistribution.allCounts(noSynapsesOnMotor,:) = NaN;
Saccadic.MotorDistributionABD = sum(Saccadic.MotorDistributionAllCounts(:,2:3),2);
Saccadic.MotorDistributionABDi = sum(Saccadic.MotorDistributionAllCounts(:,4:5),2);
Saccadic.MotorDistributionCellID = Saccadic.MotorDistributionAllCounts(:,1);
Saccadic.MotorDistributionMormalizedMotorCounts = (Saccadic.MotorDistributionABD-Saccadic.MotorDistributionABDi)./...
    (Saccadic.MotorDistributionABD+Saccadic.MotorDistributionABDi);
%SaccadicMotorDistribution.normalized(:,2) = SaccadicMotorDistribution.normalizedMotorCounts(~isnan(SaccadicMotorDistribution.normalizedMotorCounts(:,2)),2);
%SaccadicMotorDistribution.normalized(:,1) = SaccadicMotorDistribution.normalizedMotorCounts(~isnan(SaccadicMotorDistribution.normalizedMotorCounts(:,2)),1); 



% [~,sortSac] = sortrows(SaccadicMotorDistribution.normalized,2);
% cmap = colorcet('D2','N',numel(unique(SaccadicMotorDistribution.normalized(:,2))),'reverse',1);
% cmapNew = [repmat(cmap(1,:),31,1);cmap;repmat(cmap(end,:),19,1)];
% 
% transform_swc_AV(SaccadicMotorDistribution.normalized(sortSac,1),cmapNew,[],true,false);
% colormap(gca,cmapNew);
% colorbar(gca);


figure;
subplot(4,4,1)
removeZeroValues = find(Saccadic.MotorDistributionABD == 0 & Saccadic.MotorDistributionABDi ==0);
scatter(Saccadic.MotorDistributionABD(setdiff(1:length(SaccadicAxons),removeZeroValues)),...
    Saccadic.MotorDistributionABDi(setdiff(1:length(SaccadicAxons),removeZeroValues)),'ko','filled');
axis square;
set(gca,'XLim',[0,150],'ylim',[0,150]);
xlabel('Synapses on M');
ylabel('Synapses on I');
offsetAxes(gca);

subplot(4,4,2)
h = histogram2(Saccadic.MotorDistributionABD(setdiff(1:length(SaccadicAxons),removeZeroValues)),...
    Saccadic.MotorDistributionABDi(setdiff(1:length(SaccadicAxons),removeZeroValues)),...
    0:10:150,0:10:150,'DisplayStyle','tile','ShowEmptyBins','on');
axis square;
box off;
colormap(colorcet('L17','N',22,'reverse',1));
colorbar;
xlabel('Synapses onto M');
ylabel('Synapses onto I');

%SaccadicProjectingToABDexclusively = SaccadicMotorDistribution.normalizedMotorCounts(find(SaccadicMotorDistribution.normalizedMotorCounts(:,2) == 1));
%SaccadicProjectingToABDexclusivelyIndex = find(ismember(SaccadicAxons,SaccadicProjectingToABDexclusively));

SaccadicProjectingToABDexclusively.Index = find(Saccadic.MotorDistributionABD > 10 & Saccadic.MotorDistributionABDi == 0)
SaccadicProjectingToABDexclusively.cellID = SaccadicAxons(SaccadicProjectingToABDexclusively.Index);

% SaccadicProjectingToABDiexclusively = SaccadicMotorDistribution.normalizedMotorCounts(find(SaccadicMotorDistribution.normalizedMotorCounts(:,2) ==-1));
% SaccadicProjectingToABDiexclusivelyIndex = find(ismember(SaccadicAxons,SaccadicProjectingToABDiexclusively));

SaccadicProjectingToABDiexclusively.Index = find(Saccadic.MotorDistributionABD == 0 & Saccadic.MotorDistributionABDi >10);
SaccadicProjectingToABDiexclusively.cellID = SaccadicAxons(SaccadicProjectingToABDiexclusively.Index);

% SaccadicProjectingToABD_ABDi = SaccadicMotorDistribution.normalizedMotorCounts(find(SaccadicMotorDistribution.normalizedMotorCounts(:,2)>-0.95 & SaccadicMotorDistribution.normalizedMotorCounts(:,2)<0.95 ));
% SaccadicProjectingToABD_ABDi_Index = find(ismember(SaccadicAxons,SaccadicProjectingToABD_ABDi));

SaccadicProjectingToABD_ABDi.Index =  find(Saccadic.MotorDistributionABD >10 & Saccadic.MotorDistributionABDi >10);
SaccadicProjectingToABD_ABDi.cellID = SaccadicAxons(SaccadicProjectingToABD_ABDi.Index);

save('SaccadicProjectingToABDexclusively.mat','SaccadicProjectingToABDexclusively');
save('SaccadicProjectingToABDiexclusively.mat','SaccadicProjectingToABDiexclusively');
save('SaccadicProjectingToABD_ABDi.mat','SaccadicProjectingToABD_ABDi');


subplot(4,4,3)
scatter(Saccadic.MotorDistributionABD(SaccadicProjectingToABDexclusively.Index),...
    Saccadic.MotorDistributionABDi(SaccadicProjectingToABDexclusively.Index),20,SaccABDcolor,'filled');
hold on;
scatter(Saccadic.MotorDistributionABD(SaccadicProjectingToABDiexclusively.Index),...
    Saccadic.MotorDistributionABDi(SaccadicProjectingToABDiexclusively.Index),20,SaccABDicolor,'filled');
scatter(Saccadic.MotorDistributionABD(SaccadicProjectingToABD_ABDi.Index),...
    Saccadic.MotorDistributionABDi(SaccadicProjectingToABD_ABDi.Index),20,'k','filled');

axis square;
set(gca,'XLim',[0,150],'ylim',[0,150]);
xlabel('Synapses on M');
ylabel('Synapses on I');
offsetAxes(gca);

subplot(4,4,4)
clear temp;
temp = isMotor(SaccadicProjectingToABDexclusively.cellID',df);
tempSum = sum(temp(:,2:end),2);
histogram(tempSum,0:30:150,'DisplayStyle','stairs','EdgeColor',SaccABDcolor,'LineWidth',2);
hold on;
clear temp;
temp = isMotor(SaccadicProjectingToABDiexclusively.cellID',df);
tempSum = sum(temp(:,2:end),2);
histogram(tempSum,0:30:150,'DisplayStyle','stairs','EdgeColor',SaccABDicolor,'LineWidth',2);
clear temp;
temp = isMotor(SaccadicProjectingToABD_ABDi.cellID',df);
tempSum = sum(temp(:,2:end),2);
histogram(tempSum,0:30:150,'DisplayStyle','stairs','EdgeColor',[0.8,0.8,0.8],'LineWidth',2);
box off;
axis square;
legend('Sm','Si','both');
offsetAxes(gca);

subplot(4,4,5)
histogram(Saccadic.MotorDistributionMormalizedMotorCounts(SaccadicProjectingToABDexclusively.Index),...
    -1:0.1:1,'FaceColor',SaccABDcolor);
hold on;
histogram(Saccadic.MotorDistributionMormalizedMotorCounts(SaccadicProjectingToABDiexclusively.Index),...
    -1:0.1:1,'FaceColor',SaccABDicolor);
histogram(Saccadic.MotorDistributionMormalizedMotorCounts(SaccadicProjectingToABD_ABDi.Index),...
    -1:0.1:1,'FaceColor','k');
axis square;
box off;
offsetAxes(gca);


subplot(4,4,6)
histogram(Saccadic.MotorDistributionMormalizedMotorCounts(setdiff(1:length(SaccadicAxons),removeZeroValues)),...
    -1:0.1:1,'FaceColor','k');
axis square;
box off;
offsetAxes(gca);

figure('Position',[100 100 350 300]);
g = gramm('x',Saccadic.MotorDistributionABD(setdiff(1:length(SaccadicAxons),removeZeroValues)),'y',Saccadic.MotorDistributionABDi(setdiff(1:length(SaccadicAxons),removeZeroValues)));
g.geom_point();
g.geom_abline();
g.stat_cornerhist('location',150,'edges',-160:5:160,'aspect',0.5);
g.set_color_options('map',SaccABDcolor);
g.set_text_options('base_size',15);
g.set_names('x','Synapses on M','y','Synapses on I');
g.set_title('Saccadic')
g.draw();
g.export('file_name','SaccadicPLot','export_path','/Users/ashwin/Desktop','file_type','png')


figure('Position',[100 100 300 300]);
%g = gramm('x',Saccadic.MotorDistributionABD(setdiff(1:length(SaccadicAxons),removeZeroValues)),'y',Saccadic.MotorDistributionABDi(setdiff(1:length(SaccadicAxons),removeZeroValues)));
g = gramm('x',Saccadic.MotorDistributionABD(SaccadicProjectingToABDexclusively.Index),...
    'y',Saccadic.MotorDistributionABDi(SaccadicProjectingToABDexclusively.Index));
g.geom_point();
g.geom_abline();
g.stat_cornerhist('location',150,'edges',-150:10:150,'aspect',0.5);
g.set_color_options('map',SaccABDcolor);
g.set_text_options('base_size',15);
g.set_names('x','Synapses on M','y','Synapses on I');
%g.set_limit_extra( [0.05,0.05],[0.05,0.05]);
%g.axe_property('XLim',[0,250],'YLim',[0,250]);
g.draw();
g.update('x',Saccadic.MotorDistributionABD(SaccadicProjectingToABDiexclusively.Index),...
    'y',Saccadic.MotorDistributionABDi(SaccadicProjectingToABDiexclusively.Index));
g.geom_point();
%g.geom_abline();
g.stat_cornerhist('location',150,'edges',-150:10:150,'aspect',0.5);
g.set_color_options('map',SaccABDicolor);

g.set_text_options('base_size',15);
g.set_names('x','Synapses on M','y','Synapses on I');
g.draw();
g.update('x',Saccadic.MotorDistributionABD(SaccadicProjectingToABD_ABDi.Index),...
    'y',Saccadic.MotorDistributionABDi(SaccadicProjectingToABD_ABDi.Index));
g.geom_point();
%g.geom_abline();
g.stat_cornerhist('location',150,'edges',-150:10:150,'aspect',0.5);
g.set_color_options('map',[0.1,0.1,0.1]);

g.set_text_options('base_size',15);
g.set_names('x','Synapses on M','y','Synapses on I');
g.draw();







% ix1 = 1;
% ix2 = 1;
% for i = 1:size(SaccadicMotorDistribution.allCounts,1)
%     if SaccadicMotorDistribution.normalizedMotorCounts(i,2)>0.5
%     ABDSaccadicpop(ix1).cellIDs = SaccadicMotorDistribution.allCounts(i,1);
%     ABDSaccadicpop(ix1).Origin = Sacc(find(ABDSaccadicpop(ix1).cellIDs==SaccadicAxons)).Origin;
%     ABDSaccadicpop(ix1).MotorCounts = SaccadicMotorDistribution.normalizedMotorCounts(i,2);
%     ix1 = ix1+1;
%     elseif SaccadicMotorDistribution.normalizedMotorCounts(i,2)<-0.5
%      ABDiSaccadicpop(ix2).cellIDs = SaccadicMotorDistribution.allCounts(i,1);
%      ABDiSaccadicpop(ix2).Origin = Sacc(find(ABDiSaccadicpop(ix2).cellIDs==SaccadicAxons)).Origin;
%      ABDiSaccadicpop(ix2).MotorCounts = SaccadicMotorDistribution.normalizedMotorCounts(i,2);
%      ix2 = ix2+1;
%     end
% end

%  SaccadicProjectingToABD = [ABDSaccadicpop.cellIDs];
%  SaccadicProjectingToABDi = [ABDiSaccadicpop.cellIDs];
%  
 
 % make individual plots
 
%  for i = 1:size(SaccadicProjectingToABD_ABDi,1)
%      %subplot(4,4,i)
%      figure('units','normalized','outerposition',[0 0 1 1]);
%      transform_swc_AV(SaccadicProjectingToABD_ABDi(i),[0,0,0],[],true,false);
%      fname = sprintf('/Users/ashwin/Desktop/%5s.png',SaccadicProjectingToABD_ABDi(i))
%      export_fig(fname);
%      %hold on;
%      %transform_swc_AV([ABDc_CellIDs,ABDr_CellIDs],ABDcolor,[],false,false);
%      %transform_swc_AV([ABDIc_CellIDs,ABDIr_CellIDs],ABDicolor,[],false,false);
%  end


%%
SaccadicProjectingToABDexclusively.Rhombomere = isRhombomere(SaccadicProjectingToABDexclusively.cellID);
SaccadicProjectingToABDiexclusively.Rhombomere  = isRhombomere(SaccadicProjectingToABDiexclusively.cellID);
SaccadicProjectingToABD_ABDi.Rhombomere = isRhombomere(SaccadicProjectingToABD_ABDi.cellID);

% figure('Units','normalized','Position',[0,0,1,1]);
% transform_swc_AV(SaccadicProjectingToABDexclusively.cellID,SaccABDcolor,[],true,true,'r456-ABD-Vel');
% figure('Units','normalized','Position',[0,0,1,1]);
% transform_swc_AV(SaccadicProjectingToABDiexclusively.cellID,SaccABDicolor,[],true,true,'r456-ABDi-Vel');
% figure('Units','normalized','Position',[0,0,1,1]);
% transform_swc_AV(SaccadicProjectingToABD_ABDi.cellID,[0.5,0.5,0.5],[],true,true,'r456-Bino-Vel');


% figure('Units','normalized','Position',[0,0,1,1]);
% transform_swc_AV([77304,80821,80681,80743],SaccABDcolor,[],true,true,'r4-ABD-Vel');
% figure('Units','normalized','Position',[0,0,1,1]);
% transform_swc_AV([7743577437,80746,80804,80185,77684,77357],SaccABDcolor,[],true,true,'r5-ABD-Vel');
% 
% figure('Units','normalized','Position',[0,0,1,1]);
% transform_swc_AV([77645,79059,78358,76751,80285,81336,76625,81007,77434,80974,81002,80947,77656],SaccABDicolor,[],true,true,'r5-ABDi-Vel');
% figure('Units','normalized','Position',[0,0,1,1]);
% transform_swc_AV([77163,77151,77805,77460,77797,80995,77651,77132,78649,78641,81407,81027,78357,77162,80629,81395,77636,76618,77621,76627],SaccABDicolor,[],true,true,'r6-ABDi-Vel');


figure;
% Make saccadics with gradient
[Saccadic.OSIsorted,Saccadic.OSIindex] = sort(Saccadic.MotorDistributionMormalizedMotorCounts);
a = sum(Saccadic.OSIsorted == -1);
b = sum(Saccadic.OSIsorted == 1);
cmapSacc = flip(hex2rgb(SaccHEXgradient));
transform_swc_AV(SaccadicAxons(Saccadic.MotorDistributionMormalizedMotorCounts == -1)',cmapSacc(1,:),[],true,false);
transform_swc_AV(SaccadicAxons(Saccadic.MotorDistributionMormalizedMotorCounts == 1)',cmapSacc(end,:),[],false,false);
transform_swc_AV(SaccadicAxons(Saccadic.MotorDistributionMormalizedMotorCounts~=1 &...
    Saccadic.MotorDistributionMormalizedMotorCounts~=-1 & ~isnan(Saccadic.MotorDistributionMormalizedMotorCounts))',...
    cmapSacc(2:end-1,:),[],false,false);
colormap(cmapSacc);caxis([-1,1]);
colorbar;


H = hdf5read('/Users/ashwin/Documents/LabWork/SynapseDetector/contact_area.h5','/matrix');
HLables = hdf5read('/Users/ashwin/Documents/LabWork/SynapseDetector/contact_area.h5','/neuron_id_list');
Horder = ismember(HLables,SaccadicAxons);
Hsub = H(Horder(1:find(HLables==80315)),:); % till location of last saccadic axon
Saccadic.ABDcontactArea = sum(Hsub(:,find(HLables == 82140):find(HLables == 77296)),2);
Saccadic.ABDicontactArea = sum(Hsub(:,find(HLables == 78553):end),2);

NormalizedSaccadicContactArea = (Saccadic.ABDcontactArea-Saccadic.ABDicontactArea)./...
                                (Saccadic.ABDcontactArea+Saccadic.ABDicontactArea);
figure;
subplot(4,4,1)
histogram(Saccadic.MotorDistributionMormalizedMotorCounts(setdiff(1:127,removeZeroValues)),-1:0.1:1,'FaceColor','k','FaceAlpha',0.5);
hold on;
line([0.5,0.5],[0,40],'color','k','LineStyle',':');
line([-0.5,-0.5],[0,40],'color','k','LineStyle',':');
%histogram(NormalizedSaccadicContactArea,-1:0.1:1,'FaceColor','r','FaceAlpha',0.5);
xlabel('(A-Ai)/A+Ai)');
axis square;
box off;

figure; 
subplot(4,4,2)
scatter(Saccadic.ABDcontactArea,Saccadic.ABDicontactArea,25,'ko','filled');
                            
location = [vertcat(ABDSaccadicpop.Origin);vertcat(ABDiSaccadicpop.Origin)];
groups = [ones(size(vertcat(ABDSaccadicpop.Origin),1),1);2*ones(size(vertcat(ABDiSaccadicpop.Origin),1),1)];

figure;
subplot(4,4,[3,4])
scatterhist(location(:,1),location(:,2),'Group',groups,'color',[lightRed;lightBlue]);
set(gca,'YDir','reverse');


figure;
subplot(1,2,1)
[Saccadic.OSIsorted,Saccadic.OSIindex] = sort(Saccadic.MotorDistributionMormalizedMotorCounts);
a = sum(Saccadic.OSIsorted == -1);
b = sum(Saccadic.OSIsorted == 1);
cmapSacc = flip(hex2rgb(SaccHEXgradient));
transform_swc_AV(SaccadicAxons(Saccadic.MotorDistributionMormalizedMotorCounts == -1)',cmapSacc(1,:),[],true,false);
transform_swc_AV(SaccadicAxons(Saccadic.MotorDistributionMormalizedMotorCounts == 1)',cmapSacc(end,:),[],false,false);
transform_swc_AV(SaccadicAxons(Saccadic.MotorDistributionMormalizedMotorCounts~=1 &...
    Saccadic.MotorDistributionMormalizedMotorCounts~=-1 & ~isnan(Saccadic.MotorDistributionMormalizedMotorCounts))',...
    cmapSacc(2:end-1,:),[],false,false);
colormap(cmapSacc);caxis([-1,1]);
colorbar;

subplot(1,2,2)
Saccadic.Origin = getOrigin(Saccadic.MotorDistributionCellID);
scatter(Saccadic.MotorDistributionMormalizedMotorCounts,Saccadic.Origin(:,2),'o');
 %% potential Synapses analysis
AllABDPreSynapticLocaions = [vertcat(ABDr.PreSynCoordsTransformed); vertcat(ABDc.PreSynCoordsTransformed)];
AllABDiPreSynapticLoactions = [vertcat(ABDIr.PreSynCoordsTransformed); vertcat(ABDIc.PreSynCoordsTransformed)];

Saccadic.PotentialABDSynapses = zeros(numel(SaccadicAxons ),2);
Saccadic.PotentialABDiSynapses = zeros(numel(SaccadicAxons),2);


for i = 1:numel(SaccadicAxons)
    if ~isempty(Sacc(i).cellID)
        SaccABDdistanceMatrix = pdist2(Sacc(i).PostSynCoordsTransformed,AllABDPreSynapticLocaions);
        Saccadic.PotentialABDSynapses(i,1) = numel(find(SaccABDdistanceMatrix <=0.5));
        temp = isMotor(Sacc(i).cellID,df);
        Saccadic.PotentialABDSynapses(i,2) = sum(temp(2:3));
        clear temp;
        clear SaccABDdistanceMatrix;
        
        SaccABDIdistanceMatrix = pdist2(Sacc(i).PostSynCoordsTransformed,AllABDiPreSynapticLoactions);
        Saccadic.PotentialABDiSynapses(i,1) = numel(find(SaccABDIdistanceMatrix <=0.5));
        temp = isMotor(Sacc(i).cellID,df);
        Saccadic.PotentialABDiSynapses(i,2) = sum(temp(4:5));
        clear temp;
        clear Saccadic.ABDipotentialSynapses;   
    else
        Saccadic.PotentialABDSynapses(i,:) = [NaN, NaN];
        Saccadic.PotentialABDiSynapses(i,:) = [NaN, NaN];
        
    end
end

removeZeroValues = find(Saccadic.MotorDistributionABD == 0 & Saccadic.MotorDistributionABDi ==0);

subplot(4,4,1);
temp = find(ismember(SaccadicAxons,SaccadicProjectingToABDexclusively.cellID));
scatter(Saccadic.PotentialABDSynapses(temp,2),Saccadic.PotentialABDiSynapses(temp,1),30,SaccABDcolor,'filled');
hold on;
temp2 = find(ismember(SaccadicAxons,SaccadicProjectingToABDiexclusively.cellID));
scatter(Saccadic.PotentialABDiSynapses(temp2,2),Saccadic.PotentialABDSynapses(temp2,1),30,SaccABDicolor,'filled');
temp3 = find(ismember(SaccadicAxons,SaccadicProjectingToABD_ABDi.cellID));
scatter(Saccadic.PotentialABDiSynapses(temp3,2),Saccadic.PotentialABDSynapses(temp3,1),30,[0.5,0.5,0.5],'filled');

%daspect([4,1,1]);
box off;
clear temp;
clear temp2;
offsetAxes(gca);
xlabel('Actual synapses');
ylabel('Potential synapses');

%% Saccadic Motor partners

Saccadic.SortedAxons = [SaccadicProjectingToABDexclusively.cellID';
                       SaccadicProjectingToABD_ABDi.cellID';  
                       SaccadicProjectingToABDiexclusively.cellID'];


% for i = 1:numel([ABDSaccadicpop.cellIDs])
%   ABDSaccadicpop(i).MotorNeuronCounts = isPostSynapseMotor(ABDSaccadicpop(i).cellIDs,df); 
% end
% 
% for i = 1:numel([ABDiSaccadicpop.cellIDs])
%   ABDiSaccadicpop(i).MotorNeuronCounts = isPostSynapseMotor(ABDiSaccadicpop(i).cellIDs,df); 
% end

for i = 1:numel(Saccadic.SortedAxons)
    temp = isPostSynapseMotor(Saccadic.SortedAxons(i),df);
    ABDrSaccadicpopMotorConn(i,:) = temp(1:14);
    ABDcSaccadicpopMotorConn(i,:) = temp(2,1:17);
    ABDIrSaccadicpopMotorConn(i,:) = temp(3,1:11);
    ABDIcSaccadicpopMotorConn(i,:) = temp(4,1:10);
    clear temp;
end


subplot(4,4,[1,5,9])
heatmap(ABDrSaccadicpopMotorConn,'ColorScaling','scaledcolumns');

subplot(4,4,[2,6,10])
heatmap(ABDcSaccadicpopMotorConn,'ColorScaling','scaledcolumns');


subplot(4,4,[3,7,11])
heatmap(ABDIrSaccadicpopMotorConn,'ColorScaling','scaledcolumns');


subplot(4,4,[4,8,12])
heatmap(ABDIcSaccadicpopMotorConn,'ColorScaling','scaledcolumns');


subplot(4,4,13)
plot(ABDr_vols,sum(ABDrSaccadicpopMotorConn),'o','MarkerFaceColor','k','MarkerEdgeColor','none');
box off;

subplot(4,4,14)
plot(ABDc_vols,sum(ABDcSaccadicpopMotorConn),'o','MarkerFaceColor','k','MarkerEdgeColor','none');
box off;

subplot(4,4,15)
plot(ABDIr_vols,sum(ABDIrSaccadicpopMotorConn),'o','MarkerFaceColor','k','MarkerEdgeColor','none');
box off;

subplot(4,4,16)
plot(ABDIc_vols,sum(ABDIcSaccadicpopMotorConn),'o','MarkerFaceColor','k','MarkerEdgeColor','none');
box off;

%% where on the ABD are the Saccadic synapses.

load ABDr.mat
load ABDc.mat
load ABDIr.mat
load ABDIc.mat


[ABDSaccPathLengthABDr,ABDrSaccadicGradient] = getABDgradient(ABDr,vertcat(SaccadicProjectingToABDexclusively.cellID',SaccadicProjectingToABDiexclusively.cellID'),true);
[ABDSaccPathLengthABDc,ABDcSaccadicGradient] = getABDgradient(ABDc,vertcat(SaccadicProjectingToABDexclusively.cellID',SaccadicProjectingToABDiexclusively.cellID'),true);
[ABDSaccPathLengthABDIr,ABDIrSaccadicGradient] = getABDgradient(ABDIr,vertcat(SaccadicProjectingToABDexclusively.cellID',SaccadicProjectingToABDiexclusively.cellID'),true);
[ABDSaccPathLengthABDIc,ABDIcSaccadicGradient] = getABDgradient(ABDIc,vertcat(SaccadicProjectingToABDexclusively.cellID',SaccadicProjectingToABDiexclusively.cellID'),true);

[ABD_ABDi_pathLenghtABDr,ABD_ABDi_ABDrSaccadicGradient] = getABDgradient(ABDr,SaccadicProjectingToABD_ABDi.cellID ,true);
[ABD_ABDi_pathLenghtABDc,ABD_ABDi_ABDcSaccadicGradient] = getABDgradient(ABDc,SaccadicProjectingToABD_ABDi.cellID,true);
[ABD_ABDi_pathLenghtABDIr,ABD_ABDi_ABDIrSaccadicGradient] = getABDgradient(ABDIr,SaccadicProjectingToABD_ABDi.cellID,true);
[ABD_ABDi_pathLenghtABDIc,ABD_ABDi_ABDIcSaccadicGradient] = getABDgradient(ABDIc,SaccadicProjectingToABD_ABDi.cellID,true);


% figure;
% subplot(4,4,[1,5,9,13])
% imagesc(ABDrSaccadicGradient);
% set(gca,'XTickLabel',[0.2,0.4,0.6,0.8,1]);
% colormap(gca,colorcet('L8'));
% colorbar(gca);
% 
% subplot(4,4,2)
% histogram(ABDSaccPathLengthABDr{2},0:0.1:1,'Normalization','probability');
% title('row:2')
% box off;
% subplot(4,4,6)
% shadedErrorBar(0:0.1:0.9,mean(ABDrSaccadicGradient(1:7,:)),std(ABDrSaccadicGradient(1:7,:))/sqrt(7),'lineprops',{'Color','k'});
% hold on;
% shadedErrorBar(0:0.1:0.9,mean(ABDrSaccadicGradient(8:14,:)),std(ABDrSaccadicGradient(8:14,:))/sqrt(7),'lineprops',{'Color','r'});
% box off;
% subplot(4,4,14)
% histogram(ABDSaccPathLengthABDr{13},0:0.1:1,'Normalization','probability');
% title('row:13');
% box off;
% 
% subplot(4,4,[3,7,11,15])
% imagesc(ABDcSaccadicGradient);
% set(gca,'XTickLabel',[0.2,0.4,0.6,0.8,1]);
% colormap(gca,colorcet('L8'));
% colorbar(gca);
% 
% subplot(4,4,4)
% histogram(ABDSaccPathLengthABDc{2},0:0.1:1,'Normalization','probability');
% title('row:2');
% box off;
% subplot(4,4,8)
% shadedErrorBar(0:0.1:0.9,nanmean(ABDcSaccadicGradient(1:9,:)),nanstd(ABDcSaccadicGradient(1:9,:))/sqrt(8),'lineprops',{'Color','k'});
% hold on;
% shadedErrorBar(0:0.1:0.9,nanmean(ABDcSaccadicGradient(10:17,:)),nanstd(ABDcSaccadicGradient(10:17,:))/sqrt(7),'lineprops',{'Color','r'});
% box off;
% subplot(4,4,16)
% histogram(ABDSaccPathLengthABDc{15},0:0.1:1,'Normalization','probability');
% title('row:15');
% box off;
% 
% figure;
% 
% subplot(4,4,[1,5,9,13])
% imagesc(ABDIrSaccadicGradient);
% set(gca,'XTickLabel',[0.2,0.4,0.6,0.8,1]);
% colormap(gca,colorcet('L8'));
% colorbar(gca);
% subplot(4,4,2)
% histogram(ABDSaccPathLengthABDIr{3},0:0.1:1,'Normalization','probability');
% title('row:3');
% box off;
% subplot(4,4,6)
% shadedErrorBar(0:0.1:0.9,mean(ABDIrSaccadicGradient(1:5,:)),std(ABDIrSaccadicGradient(1:5,:))/sqrt(5),'lineprops',{'Color','k'});
% hold on;
% shadedErrorBar(0:0.1:0.9,mean(ABDIrSaccadicGradient(6:11,:)),std(ABDIrSaccadicGradient(6:11,:))/sqrt(6),'lineprops',{'Color','r'});
% box off;
% subplot(4,4,14)
% histogram(ABDSaccPathLengthABDIr{10},0:0.1:1,'Normalization','probability');
% title('row:10');
% box off;
% 
% subplot(4,4,[3,7,11,15])
% imagesc(ABDIcSaccadicGradient);
% set(gca,'XTickLabel',[0.2,0.4,0.6,0.8,1]);
% colormap(gca,colorcet('L8'));
% colorbar(gca);
% 
% subplot(4,4,4)
% histogram(ABDSaccPathLengthABDIc{3},0:0.1:1,'Normalization','probability');
% title('row:3');
% box off;
% subplot(4,4,8)
% shadedErrorBar(0:0.1:0.9,nanmean(ABDIcSaccadicGradient(1:5,:)),nanstd(ABDIcSaccadicGradient(1:5,:))/sqrt(5),'lineprops',{'Color','k'});
% hold on;
% shadedErrorBar(0:0.1:0.9,nanmean(ABDIcSaccadicGradient(6:10,:)),nanstd(ABDIcSaccadicGradient(6:10,:))/sqrt(4),'lineprops',{'Color','r'});
% box off;
% subplot(4,4,16)
% histogram(ABDSaccPathLengthABDIc{10},0:0.1:1,'Normalization','probability');
% title('row:10');
% box off;

%%


figure;
subplot(4,4,1)
histogram([cell2mat(ABDSaccPathLengthABDr),cell2mat(ABDSaccPathLengthABDc)],'EdgeColor','r','DisplayStyle','stairs','LineWidth',2);
hold on;
histogram([cell2mat(ABD_ABDi_pathLenghtABDr),cell2mat(ABD_ABDi_pathLenghtABDc)],'EdgeColor','b','DisplayStyle','stairs','LineWidth',2);
axis square;
box off;
xlabel('Norm. pathlength');
ylabel('counts');
title('ABD/ABD_ABDi');

subplot(4,4,2)
histogram([cell2mat(ABDSaccPathLengthABDIr),cell2mat(ABDSaccPathLengthABDIc)],'EdgeColor','r','DisplayStyle','stairs','LineWidth',2);
hold on;
histogram([cell2mat(ABD_ABDi_pathLenghtABDIr),cell2mat(ABD_ABDi_pathLenghtABDIc)],'EdgeColor','b','DisplayStyle','stairs','LineWidth',2);
axis square;
box off;
xlabel('Norm. pathlength');
ylabel('counts');
title('ABDi/ABD_ABDi');


Saccadic.meanABDgradient = nanmean(vertcat(ABDrSaccadicGradient,ABDcSaccadicGradient));
Saccadic.stdABDgradient = nanstd(vertcat(ABDrSaccadicGradient,ABDcSaccadicGradient));
Saccadic.numberOfABDneurons = 29;

Saccadic.meanABDigradient = nanmean(vertcat(ABDIrSaccadicGradient,ABDIcSaccadicGradient));
Saccadic.stdABDigradient = nanstd(vertcat(ABDIrSaccadicGradient,ABDIcSaccadicGradient));
Saccadic.numberOfABDineurons = 21;

Saccadic.meanABD_ABDi_ABDgradient = nanmean(vertcat(ABD_ABDi_ABDrSaccadicGradient,ABD_ABDi_ABDcSaccadicGradient));
Saccadic.stdABD_ABDi_ABDgradient = nanstd(vertcat(ABD_ABDi_ABDrSaccadicGradient,ABD_ABDi_ABDcSaccadicGradient));

Saccadic.meanABD_ABDi_ABDigradient = nanmean(vertcat(ABD_ABDi_ABDIrSaccadicGradient,ABD_ABDi_ABDIcSaccadicGradient));
Saccadic.stdABD_ABDi_ABDigradient = nanstd(vertcat(ABD_ABDi_ABDIrSaccadicGradient,ABD_ABDi_ABDIcSaccadicGradient));


figure('Units','normalized','Position',[0,0,1,1]);
subplot(4,4,1)
errorbar(Saccadic.meanABDgradient,Saccadic.stdABDgradient./sqrt(Saccadic.numberOfABDneurons),...
    '-o','color',SaccABDcolor,'LineWidth',2,'MarkerFaceColor','w');
hold on;
errorbar(Saccadic.meanABD_ABDi_ABDgradient,Saccadic.stdABD_ABDi_ABDgradient./sqrt(Saccadic.numberOfABDneurons),...
    '-o','color',[0.5,0.5,0.5],'LineWidth',2,'MarkerFaceColor','w');
axis square;
set(gca,'XTickLabels',[0,0.5,1],'color',[ABDcolor,0.1]);
xlabel('Norm. pathlength');
ylabel('Norm. count');
box off;
offsetAxes(gca);

subplot(4,4,2)
errorbar(Saccadic.meanABDigradient,Saccadic.stdABDigradient./sqrt(Saccadic.numberOfABDineurons),...
    '-o','color',SaccABDicolor,'LineWidth',2,'MarkerFaceColor','w');
hold on;
errorbar(Saccadic.meanABD_ABDi_ABDigradient,Saccadic.stdABD_ABDi_ABDigradient./sqrt(Saccadic.numberOfABDineurons),...
    '-o','color',[0.5,0.5,0.5],'LineWidth',2,'MarkerFaceColor','w');

axis square;
set(gca,'XTickLabels',[0,0.5,1],'color',[ABDicolor,0.1]);
xlabel('Norm. pathlength');
ylabel('Norm. count');
box off;
offsetAxes(gca);


subplot(4,4,5)
histogram(vertcat(ABDrSaccadicGradient,ABDcSaccadicGradient),20,'Normalization','cdf',...
    'displayStyle','stairs','EdgeColor',SaccABDcolor,'LineWidth',2);
hold on;
histogram(vertcat(ABD_ABDi_ABDrSaccadicGradient,ABD_ABDi_ABDcSaccadicGradient),20,'Normalization','cdf',...
    'displayStyle','stairs','EdgeColor',[0.5,0.5,0.5],'LineWidth',2);
box off;
axis square;
xlabel('Norm. pathlength');
ylabel('Cumulative count');
offsetAxes(gca);

subplot(4,4,6)
histogram(vertcat(ABDIrSaccadicGradient,ABDIcSaccadicGradient),20,'Normalization','cdf',...
    'displayStyle','stairs','EdgeColor',SaccABDicolor,'LineWidth',2);
hold on;
histogram(vertcat(ABD_ABDi_ABDIrSaccadicGradient,ABD_ABDi_ABDIcSaccadicGradient),20,'Normalization','cdf',...
    'displayStyle','stairs','EdgeColor',[0.5,0.5,0.5],'LineWidth',2);
box off;
axis square;
xlabel('Norm. pathlength');
ylabel('Cumulative count');
offsetAxes(gca);



%% Sort Saccadics based on rhombomere position.

SaccadicProjectingToABDexclusively.origin = vertcat(Sacc(SaccadicProjectingToABDexclusively.Index).Origin);
[~,SaccadicProjectingToABDexclusively.RCsort] = sortrows(SaccadicProjectingToABDexclusively.origin,2);
%SaccadicProjectingToABDexclusively.r4 = SaccadicProjectingToABDexclusively.cellID(SaccadicProjectingToABDexclusively.origin(:,2)<625);
%SaccadicProjectingToABDexclusively.r5 = SaccadicProjectingToABDexclusively.cellID(SaccadicProjectingToABDexclusively.origin(:,2)>625 & SaccadicProjectingToABDexclusively.origin(:,2)<660);
SaccadicProjectingToABDexclusively.r6 = SaccadicProjectingToABDexclusively.cellID(SaccadicProjectingToABDexclusively.origin(:,2)>660);
SaccadicProjectingToABDexclusively.r4 = [80821,80217,80801,80763,77373,80743,78346,80850,77342,77304,80681];
SaccadicProjectingToABDexclusively.r5 = [77374,77357,77435,80804,77461,77684,77433,80746,77437,80185];

SaccadicProjectingToABDiexclusively.origin = vertcat(Sacc(SaccadicProjectingToABDiexclusively.Index).Origin);
[~,SaccadicProjectingToABDiexclusively.RCsort] = sortrows(SaccadicProjectingToABDiexclusively.origin,2);
SaccadicProjectingToABDiexclusively.r4 = SaccadicProjectingToABDiexclusively.cellID(SaccadicProjectingToABDiexclusively.origin(:,2)<625);
%SaccadicProjectingToABDiexclusively.r5 = SaccadicProjectingToABDiexclusively.cellID(SaccadicProjectingToABDiexclusively.origin(:,2)>625 & SaccadicProjectingToABDiexclusively.origin(:,2)<660);
%SaccadicProjectingToABDiexclusively.r6 = SaccadicProjectingToABDiexclusively.cellID(SaccadicProjectingToABDiexclusively.origin(:,2)>660);
SaccadicProjectingToABDiexclusively.r5  = [81406,81007,76622,77645,79059,78650,81002,77434,78358,76625,80974,76751,77656,78544];
SaccadicProjectingToABDiexclusively.r6  = [78641,77805,77122,77162,77132,78357,77797,77460,76627,77651,77151,77163,77636,76697,76618,77848,80995,77390,77621];

SaccadicProjectingToABD_ABDi.rhombomere = vertcat(Sacc(SaccadicProjectingToABD_ABDi.Index).Rhombomere);
SaccadicProjectingToABD_ABDi.origin = vertcat(Sacc(SaccadicProjectingToABD_ABDi.Index).Origin);
[~,SaccadicProjectingToABD_ABDi.RCsort] = sortrows(SaccadicProjectingToABD_ABDi.origin,2);


SaccadicProjectingToABDexclusively.motorDistributionAllCounts = Saccadic.MotorDistributionAllCounts(SaccadicProjectingToABDexclusively.Index,:);
SaccadicProjectingToABDexclusively.motorDistributionABD = Saccadic.MotorDistributionABD(SaccadicProjectingToABDexclusively.Index);
SaccadicProjectingToABDexclusively.motorDistributionABDi = Saccadic.MotorDistributionABDi(SaccadicProjectingToABDexclusively.Index);

SaccadicProjectingToABDiexclusively.motorDistributionAllCounts = Saccadic.MotorDistributionAllCounts(SaccadicProjectingToABDiexclusively.Index,:);
SaccadicProjectingToABDiexclusively.motorDistributionABD = Saccadic.MotorDistributionABD(SaccadicProjectingToABDiexclusively.Index);
SaccadicProjectingToABDiexclusively.motorDistributionABDi = Saccadic.MotorDistributionABDi(SaccadicProjectingToABDiexclusively.Index);

SaccadicProjectingToABD_ABDi.motorDistributionAllCounts = Saccadic.MotorDistributionAllCounts(SaccadicProjectingToABD_ABDi.Index,:);
SaccadicProjectingToABD_ABDi.motorDistributionABD = Saccadic.MotorDistributionABD(SaccadicProjectingToABD_ABDi.Index);
SaccadicProjectingToABD_ABDi.motorDistributionABDi = Saccadic.MotorDistributionABDi(SaccadicProjectingToABD_ABDi.Index);


figure;

subplot(3,3,[1,4])
heatmap(SaccadicProjectingToABDexclusively.motorDistributionAllCounts(SaccadicProjectingToABDexclusively.RCsort,2:5),...
    'ColorScaling','scaledcolumns');

subplot(3,3,[2,5])
heatmap(SaccadicProjectingToABDiexclusively.motorDistributionAllCounts(SaccadicProjectingToABDiexclusively.RCsort,2:5),...
    'ColorScaling','scaledcolumns');

subplot(3,3,3)
heatmap(SaccadicProjectingToABD_ABDi.motorDistributionAllCounts(SaccadicProjectingToABD_ABDi.RCsort,2:5),...
    'ColorScaling','scaledcolumns');

subplot(3,3,7)
scatter(SaccadicProjectingToABDexclusively.motorDistributionAllCounts(:,2),SaccadicProjectingToABDexclusively.motorDistributionAllCounts(:,3),'o',...
   'MarkerFaceColor', SaccABDcolor,'MarkerEdgeColor','none');
hold on;
%showfit(ezfit(SaccadicProjectingToABDexclusively.motorDistributionAllCounts(:,2),SaccadicProjectingToABDexclusively.motorDistributionAllCounts(:,3),'a*x+b'),'fitcolor',SaccABDcolor,'dispeqboxmode','off','dispfitlegend','on','corrcoefmode','r2');
scatter(SaccadicProjectingToABDiexclusively.motorDistributionAllCounts(:,4),SaccadicProjectingToABDiexclusively.motorDistributionAllCounts(:,5),'o',...
   'MarkerFaceColor', SaccABDicolor,'MarkerEdgeColor','none');
%showfit(ezfit(SaccadicProjectingToABDiexclusively.motorDistributionAllCounts(:,4),SaccadicProjectingToABDiexclusively.motorDistributionAllCounts(:,5),'a*x+b'),'fitcolor',SaccABDicolor,'dispeqboxmode','off','dispfitlegend','on','corrcoefmode','r2');
xlabel('Synapses on ABDr/ABDIr');
ylabel('Synapses on ABDc/ABDIc');
axis square;



% motor gradient patterns

% r4

for i = 1:numel(ABDr_CellIDs)
    temp1 = ismember(ABDr(i).Inputs,SaccadicProjectingToABDexclusively.r4');
    temp2 = ABDr(i).PathLength(temp1)'./max(Pvec_tree(ABDr(i).Tree{1}));
    SaccadicProjectingToABDexclusively.ABDrR4Gradient(i,:) =  histcounts(temp2,0:0.1:1);
    clear temp1;
    clear temp2;
end

for i = 1:numel(ABDc_CellIDs)
    if ~isempty(ABDc(i).Tree)
        temp1 = ismember(ABDc(i).Inputs,SaccadicProjectingToABDexclusively.r4');
        temp2 = ABDc(i).PathLength(temp1)'./max(Pvec_tree(ABDc(i).Tree{1}));
        SaccadicProjectingToABDexclusively.ABDcR4Gradient(i,:) =  histcounts(temp2,0:0.1:1);
        clear temp1;
        clear temp2;
    end
end

for i = 1:numel(ABDIr_CellIDs)
    temp1 = ismember(ABDIr(i).Inputs,SaccadicProjectingToABDexclusively.r4');
    temp2 = ABDIr(i).PathLength(temp1)'./max(Pvec_tree(ABDIr(i).Tree{1}));
    SaccadicProjectingToABDexclusively.ABDIrR4Gradient(i,:) =  histcounts(temp2,0:0.1:1);
    clear temp1;
    clear temp2;
end

for i = 1:numel(ABDIc_CellIDs)
    if ~isempty(ABDIc(i).Tree)
        temp1 = ismember(ABDIc(i).Inputs,SaccadicProjectingToABDexclusively.r4');
        temp2 = ABDIc(i).PathLength(temp1)'./max(Pvec_tree(ABDIc(i).Tree{1}));
        SaccadicProjectingToABDexclusively.ABDIcR4Gradient(i,:) =  histcounts(temp2,0:0.1:1);
        clear temp1;
        clear temp2;
    end
end

%r5


for i = 1:numel(ABDr_CellIDs)
    temp1 = ismember(ABDr(i).Inputs,SaccadicProjectingToABDexclusively.r5');
    temp2 = ABDr(i).PathLength(temp1)'./max(Pvec_tree(ABDr(i).Tree{1}));
    SaccadicProjectingToABDexclusively.ABDrR5Gradient(i,:) =  histcounts(temp2,0:0.1:1);
    clear temp1;
    clear temp2;
end

for i = 1:numel(ABDc_CellIDs)
    if ~isempty(ABDc(i).Tree)
        temp1 = ismember(ABDc(i).Inputs,SaccadicProjectingToABDexclusively.r5');
        temp2 = ABDc(i).PathLength(temp1)'./max(Pvec_tree(ABDc(i).Tree{1}));
        SaccadicProjectingToABDexclusively.ABDcR5Gradient(i,:) =  histcounts(temp2,0:0.1:1);
        clear temp1;
        clear temp2;
    end
end

for i = 1:numel(ABDIr_CellIDs)
    temp1 = ismember(ABDIr(i).Inputs,SaccadicProjectingToABDexclusively.r5');
    temp2 = ABDIr(i).PathLength(temp1)'./max(Pvec_tree(ABDIr(i).Tree{1}));
    SaccadicProjectingToABDexclusively.ABDIrR5Gradient(i,:) =  histcounts(temp2,0:0.1:1);
    clear temp1;
    clear temp2;
end

for i = 1:numel(ABDIc_CellIDs)
    if ~isempty(ABDIc(i).Tree)
        temp1 = ismember(ABDIc(i).Inputs,SaccadicProjectingToABDexclusively.r5');
        temp2 = ABDIc(i).PathLength(temp1)'./max(Pvec_tree(ABDIc(i).Tree{1}));
        SaccadicProjectingToABDexclusively.ABDIcR5Gradient(i,:) =  histcounts(temp2,0:0.1:1);
        clear temp1;
        clear temp2;
    end
end

%r6

for i = 1:numel(ABDr_CellIDs)
    temp1 = ismember(ABDr(i).Inputs,SaccadicProjectingToABDexclusively.r6');
    temp2 = ABDr(i).PathLength(temp1)'./max(Pvec_tree(ABDr(i).Tree{1}));
    SaccadicProjectingToABDexclusively.ABDrR6Gradient(i,:) =  histcounts(temp2,0:0.1:1);
    clear temp1;
    clear temp2;
end

for i = 1:numel(ABDc_CellIDs)
    if ~isempty(ABDc(i).Tree)
        temp1 = ismember(ABDc(i).Inputs,SaccadicProjectingToABDexclusively.r6');
        temp2 = ABDc(i).PathLength(temp1)'./max(Pvec_tree(ABDc(i).Tree{1}));
        SaccadicProjectingToABDexclusively.ABDcR6Gradient(i,:) =  histcounts(temp2,0:0.1:1);
        clear temp1;
        clear temp2;
    end
end

for i = 1:numel(ABDIr_CellIDs)
    temp1 = ismember(ABDIr(i).Inputs,SaccadicProjectingToABDexclusively.r6');
    temp2 = ABDIr(i).PathLength(temp1)'./max(Pvec_tree(ABDIr(i).Tree{1}));
    SaccadicProjectingToABDexclusively.ABDIrR6Gradient(i,:) =  histcounts(temp2,0:0.1:1);
    clear temp1;
    clear temp2;
end

for i = 1:numel(ABDIc_CellIDs)
    if ~isempty(ABDc(i).Tree)
        temp1 = ismember(ABDIc(i).Inputs,SaccadicProjectingToABDexclusively.r6');
        temp2 = ABDIc(i).PathLength(temp1)'./max(Pvec_tree(ABDIc(i).Tree{1}));
        SaccadicProjectingToABDexclusively.ABDIcR6Gradient(i,:) =  histcounts(temp2,0:0.1:1);
        clear temp1;
        clear temp2;
    end
end

%ABDi

%r5

for i = 1:numel(ABDr_CellIDs)
    temp1 = ismember(ABDr(i).Inputs,SaccadicProjectingToABDiexclusively.r5');
    temp2 = ABDr(i).PathLength(temp1)'./max(Pvec_tree(ABDr(i).Tree{1}));
    SaccadicProjectingToABDiexclusively.ABDrR5Gradient(i,:) =  histcounts(temp2,0:0.1:1);
    clear temp1;
    clear temp2;
end

for i = 1:numel(ABDc_CellIDs)
    if ~isempty(ABDc(i).Tree)
        temp1 = ismember(ABDc(i).Inputs,SaccadicProjectingToABDiexclusively.r5');
        temp2 = ABDc(i).PathLength(temp1)'./max(Pvec_tree(ABDc(i).Tree{1}));
        SaccadicProjectingToABDiexclusively.ABDcR5Gradient(i,:) =  histcounts(temp2,0:0.1:1);
        clear temp1;
        clear temp2;
    end
end

for i = 1:numel(ABDIr_CellIDs)
    temp1 = ismember(ABDIr(i).Inputs,SaccadicProjectingToABDiexclusively.r5');
    temp2 = ABDIr(i).PathLength(temp1)'./max(Pvec_tree(ABDIr(i).Tree{1}));
    SaccadicProjectingToABDiexclusively.ABDIrR5Gradient(i,:) =  histcounts(temp2,0:0.1:1);
    clear temp1;
    clear temp2;
end

for i = 1:numel(ABDIc_CellIDs)
    if ~isempty(ABDIc(i).Tree)
        temp1 = ismember(ABDIc(i).Inputs,SaccadicProjectingToABDiexclusively.r5');
        temp2 = ABDIc(i).PathLength(temp1)'./max(Pvec_tree(ABDIc(i).Tree{1}));
        SaccadicProjectingToABDiexclusively.ABDIcR5Gradient(i,:) =  histcounts(temp2,0:0.1:1);
        clear temp1;
        clear temp2;
    end
end

%r6


for i = 1:numel(ABDr_CellIDs)
    temp1 = ismember(ABDr(i).Inputs,SaccadicProjectingToABDiexclusively.r6');
    temp2 = ABDr(i).PathLength(temp1)'./max(Pvec_tree(ABDr(i).Tree{1}));
    SaccadicProjectingToABDiexclusively.ABDrR6Gradient(i,:) =  histcounts(temp2,0:0.1:1);
    clear temp1;
    clear temp2;
end

for i = 1:numel(ABDc_CellIDs)
    if ~isempty(ABDc(i).Tree)
        temp1 = ismember(ABDc(i).Inputs,SaccadicProjectingToABDiexclusively.r6');
        temp2 = ABDc(i).PathLength(temp1)'./max(Pvec_tree(ABDc(i).Tree{1}));
        SaccadicProjectingToABDiexclusively.ABDcR6Gradient(i,:) =  histcounts(temp2,0:0.1:1);
        clear temp1;
        clear temp2;
    end
end

for i = 1:numel(ABDIr_CellIDs)
    temp1 = ismember(ABDIr(i).Inputs,SaccadicProjectingToABDiexclusively.r6');
    temp2 = ABDIr(i).PathLength(temp1)'./max(Pvec_tree(ABDIr(i).Tree{1}));
    SaccadicProjectingToABDiexclusively.ABDIrR6Gradient(i,:) =  histcounts(temp2,0:0.1:1);
    clear temp1;
    clear temp2;
end

for i = 1:numel(ABDIc_CellIDs)
    if ~isempty(ABDIc(i).Tree)
        temp1 = ismember(ABDIc(i).Inputs,SaccadicProjectingToABDiexclusively.r6');
        temp2 = ABDIc(i).PathLength(temp1)'./max(Pvec_tree(ABDIc(i).Tree{1}));
        SaccadicProjectingToABDiexclusively.ABDIcR6Gradient(i,:) =  histcounts(temp2,0:0.1:1);
        clear temp1;
        clear temp2;
    end
end

figure;

subplot(3,3,[1,4,7])
scatter3(SaccadicProjectingToABDexclusively.origin(:,1),SaccadicProjectingToABDexclusively.origin(:,2),SaccadicProjectingToABDexclusively.origin(:,3),...
    'markerFaceColor',SaccABDcolor);
hold on;
scatter3(SaccadicProjectingToABDiexclusively.origin(:,1),SaccadicProjectingToABDiexclusively.origin(:,2),SaccadicProjectingToABDiexclusively.origin(:,3),...
    'markerFaceColor',SaccABDicolor);
set(gca,'Ydir','reverse');

line([250,290],[625,625],'color','k');
line([250,290],[660,660],'color','k');
xlabel('M <--> L');
ylabel('C <--> R ');

legend('Sac-->ABD','Sac-->ABDi');
title('Somata pos');
daspect([1,1,1]);
view([0,90]);

subplot(3,3,2)
%r4
plot(sum(vertcat(SaccadicProjectingToABDexclusively.ABDrR4Gradient,SaccadicProjectingToABDexclusively.ABDcR4Gradient)),...
    '-o','MarkerFaceColor',SaccABDcolor);
title('Synapses onto ABD');

subplot(3,3,3)
%r4
plot(sum(vertcat(SaccadicProjectingToABDexclusively.ABDIrR4Gradient,SaccadicProjectingToABDexclusively.ABDIcR4Gradient)),...
    '-o','MarkerFaceColor',SaccABDcolor);
title('Synapses onto ABDi');


subplot(3,3,5)

plot(sum(vertcat(SaccadicProjectingToABDexclusively.ABDrR5Gradient,SaccadicProjectingToABDexclusively.ABDcR5Gradient)),...
    '-o','MarkerFaceColor',SaccABDcolor);
hold on;
plot(sum(vertcat(SaccadicProjectingToABDiexclusively.ABDrR5Gradient,SaccadicProjectingToABDiexclusively.ABDcR5Gradient)),...
    '-o','MarkerFaceColor',SaccABDicolor);

subplot(3,3,6)

plot(sum(vertcat(SaccadicProjectingToABDexclusively.ABDIrR5Gradient,SaccadicProjectingToABDexclusively.ABDIcR5Gradient)),...
    '-o','MarkerFaceColor',SaccABDcolor);
hold on;
plot(sum(vertcat(SaccadicProjectingToABDiexclusively.ABDIrR5Gradient,SaccadicProjectingToABDiexclusively.ABDIcR5Gradient)),...
    '-o','MarkerFaceColor',SaccABDicolor);



subplot(3,3,8)

plot(sum(vertcat(SaccadicProjectingToABDexclusively.ABDrR6Gradient,SaccadicProjectingToABDexclusively.ABDcR6Gradient)),...
    '-o','MarkerFaceColor',SaccABDcolor);
hold on;
plot(sum(vertcat(SaccadicProjectingToABDiexclusively.ABDrR6Gradient,SaccadicProjectingToABDiexclusively.ABDcR6Gradient)),...
    '-o','MarkerFaceColor',SaccABDicolor);
xlabel('ABD Norm. pathlength');
ylabel('No. synapses');

subplot(3,3,9)

plot(sum(vertcat(SaccadicProjectingToABDexclusively.ABDIrR6Gradient,SaccadicProjectingToABDexclusively.ABDIcR6Gradient)),...
    '-o','MarkerFaceColor',SaccABDcolor);
hold on;
plot(sum(vertcat(SaccadicProjectingToABDiexclusively.ABDIrR6Gradient,SaccadicProjectingToABDiexclusively.ABDIcR6Gradient)),...
    '-o','MarkerFaceColor',SaccABDicolor);

xlabel('ABDi Norm. pathlength');
ylabel('No. synapses');



%% Saccadic --> Integrator (without function)


ALX.motorDistribution  = isMotor(ALX.cellIDs',df);
ALX.OSI = (sum(ALX.motorDistribution(:,2:3),2) - sum(ALX.motorDistribution(:,4:5),2)) ./ ...
            (sum(ALX.motorDistribution(:,2:3),2) + sum(ALX.motorDistribution(:,4:5),2));



SaccadicToIntegrator.IntegratorCellIDs = unique(vertcat(Sacc.Integrator));
   
for i = 1:size(SaccadicAxons,2)
    SaccadicToIntegrator.SaccadicCellID(i) = SaccadicAxons(i); % Saccadic neurons that project to Integrators
    SaccadicToIntegrator.IntegratorCount(i,:) = histcounts(categorical(Sacc(i).Integrator),categorical(SaccadicToIntegrator.IntegratorCellIDs));
end

SaccadicToIntegrator.MotorDistribution = isMotor(SaccadicToIntegrator.SaccadicCellID',df);


% get counts onto integrators from axon projecting to ABD and ABDi
SaccadicToIntegrator.ABDIntCount = sum(SaccadicToIntegrator.IntegratorCount(SaccadicProjectingToABDexclusively.Index,:),1);
SaccadicToIntegrator.ABDiIntCount = sum(SaccadicToIntegrator.IntegratorCount(SaccadicProjectingToABDiexclusively.Index,:),1);
SaccadicToIntegrator.ABDIntNormalizedCount = (SaccadicToIntegrator.ABDIntCount-SaccadicToIntegrator.ABDiIntCount) ./ (SaccadicToIntegrator.ABDIntCount+SaccadicToIntegrator.ABDiIntCount);

% get integrator Origins

for i = 1:numel(SaccadicToIntegrator.IntegratorCellIDs)
    SaccadicToIntegrator.IntegratorOrigin(i,:) = getOrigin(SaccadicToIntegrator.IntegratorCellIDs(i));
end
    


figure;
subplot(4,4,1)
histogram(SaccadicToIntegrator.ABDIntNormalizedCount,-1:0.1:1,'FaceColor','k');
box off;
xlabel('(Sa-Sai)/(Sa+Sai)');
axis square;


subplot(4,4,3)
histogram(SaccadicToIntegrator.ABDIntNormalizedCount(ismember(SaccadicToIntegrator.IntegratorCellIDs,allDBX)),-1:0.1:1,'FaceColor',DBXcolor)
box off;
xlabel('(Pm-Pi)/(Pm+Pi)');
title('Dbx')
axis square;
offsetAxes(gca);

figure;
subplot(1,2,1)
SaccadicToIntegrator.SaccABDIntOrder = find(SaccadicToIntegrator.ABDIntNormalizedCount>0);
SaccadicToIntegrator.SaccABDIntCellIDs = SaccadicToIntegrator.IntegratorCellIDs(SaccadicToIntegrator.SaccABDIntOrder);
transform_swc_AV(SaccadicToIntegrator.SaccABDIntCellIDs,ALXcolor,[],true,false);

subplot(1,2,2)
SaccadicToIntegrator.SaccABDiIntOrder = find(SaccadicToIntegrator.ABDIntNormalizedCount<0);
SaccadicToIntegrator.SaccABDiIntCellIDs = SaccadicToIntegrator.IntegratorCellIDs(SaccadicToIntegrator.SaccABDiIntOrder);
transform_swc_AV(SaccadicToIntegrator.SaccABDiIntCellIDs,ALXcolor,[],true,false);

figure;
transform_swc_AV(allALX,ALXcolor,[],true,false);
figure;
transform_swc_AV(allDBX,DBXcolor,[],true,false);


figure;
imagesc(SaccadicToIntegrator.IntegratorCount([SaccadicProjectingToABDexclusively.Index;SaccadicProjectingToABDiexclusively.Index],...
    [SaccadicToIntegrator.SaccABDIntOrder';SaccadicToIntegrator.SaccABDiIntOrder']))
colormap(colorcet('L8'));
line([0,size(SaccadicToIntegrator.IntegratorCellIDs,1)],[size(SaccadicProjectingToABDexclusively.Index,1),size(SaccadicProjectingToABDexclusively.Index,1)],'color','k');
line([size(SaccadicToIntegrator.SaccABDIntOrder,2),size(SaccadicToIntegrator.SaccABDIntOrder,2)],[0,size(SaccadicToIntegrator.SaccadicCellID,2)],'color','k');

%% Do integrators obey ABD order or Sac-->ABD order?

allIntegrators = [confirmedALX,putativeALX,confirmedDBX,putativeDBX,confirmedBARHL,putativeBARHL];
Integrators.motorDistribution = isMotor(allIntegrators',df);
Integrators.NormalizedCount = (sum(Integrators.motorDistribution(:,2:3),2)-sum(Integrators.motorDistribution(:,4:5),2))./...
                               (sum(Integrators.motorDistribution(:,2:3),2)+sum(Integrators.motorDistribution(:,4:5),2)) ;
% Integrators that receive Saccadic Inputs
[~,allIntegratorOrder] = ismember(SaccadicToIntegrator.IntegratorCellIDs(ismember(SaccadicToIntegrator.IntegratorCellIDs,ALX.cellIDs)),allIntegrators);

figure;
subplot(4,4,2)
histogram(ALX.OSI(ismember(ALX.cellIDs,[SaccadicToIntegrator.SaccABDiIntCellIDs;SaccadicToIntegrator.SaccABDIntCellIDs])),...
    -1:0.1:1,'FaceColor',ALXcolor);
axis square;
box off;
set(gca,'XColor','none');
 ylabel('Count')
offsetAxes(gca);

%jitterAmount = 0.2;
subplot(4,4,6);
scatter(Integrators.NormalizedCount(allIntegratorOrder),SaccadicToIntegrator.ABDIntNormalizedCount(ismember(SaccadicToIntegrator.IntegratorCellIDs,ALX.cellIDs)),20,ALXcolor,'filled');
%set(gca,'YAxisLocation','right');
xlabel('OSI');
%ylabel('\textsf{$$(P_m-P_i)/(P_m+P_i)$$}','Interpreter','latex');
offsetAxes(gca);
axis square;
box off;

subplot(4,4,5)
histogram(SaccadicToIntegrator.ABDIntNormalizedCount(ismember(SaccadicToIntegrator.IntegratorCellIDs,ALX.motorDistribution(:,1))),-1:0.1:1,'FaceColor',ALXcolor)
box off;
 ylabel('Count');
% title('Alx');
axis square;
offsetAxes(gca);
set(gca, 'Xcolor','none');
camroll(90);





%% Saccade --> Saccade interactions
load AllCells.mat
load ConnMatrixPre.mat

% unbiased saccadic counts
for i = 1:size(SaccadicAxons,2)
    Saccadic.SaccadicCount(i,:) = histcounts(categorical(Sacc(i).Saccadic),categorical(SaccadicAxons)); % down stream neurons
    Saccadic.SaccadicCountNorm(i,:) = Saccadic.SaccadicCount(i,:)./ size(Sacc(i).Outputs,1);
end

[~,~,LocateSaccadicinAllCells] = intersect(SaccadicAxons,AllCells,'stable');
SaccadicProjectingToABDexclusively.RCsortinAllCells = LocateSaccadicinAllCells((SaccadicProjectingToABDexclusively.Index(SaccadicProjectingToABDexclusively.RCsort)));
SaccadicProjectingToABDiexclusively.RCsortinAllCells = LocateSaccadicinAllCells((SaccadicProjectingToABDiexclusively.Index(SaccadicProjectingToABDiexclusively.RCsort)));

figure;
subplot(4,4,1) ; % ordered by S->M and S->I
cspy(ConnMatrixPre([SaccadicProjectingToABDexclusively.RCsortinAllCells;SaccadicProjectingToABDiexclusively.RCsortinAllCells],...
    [SaccadicProjectingToABDexclusively.RCsortinAllCells;SaccadicProjectingToABDiexclusively.RCsortinAllCells]),...
    'Colormap',colorcet('R3','N',15),'Levels',15,'MarkerSize',7);

% imagesc(ConnMatrixPre([SaccadicProjectingToABDexclusively.RCsortinAllCells;SaccadicProjectingToABDiexclusively.RCsortinAllCells],...
%     [SaccadicProjectingToABDexclusively.RCsortinAllCells;SaccadicProjectingToABDiexclusively.RCsortinAllCells]));
axis square;

a = numel(SaccadicProjectingToABDexclusively.r4);
line([a,a]+0.5,[0,59]+0.5,'color','k');
b= numel(SaccadicProjectingToABDexclusively.r5);
line([a+b,a+b]+0.5,[0,59]+0.5,'color','k');
c = numel(SaccadicProjectingToABDexclusively.r6);
line([a+b+c,a+b+c]+0.5,[0,59]+0.5,'color','k');
d = numel(SaccadicProjectingToABDiexclusively.r5);
line([a+b+c+d,a+b+c+d]+0.5,[0,59]+0.5,'color','k');

line([0,59]+0.5,[a,a]+0.5,'color','k');
line([0,59]+0.5,[a+b,a+b]+0.5,'color','k');
line([0,59]+0.5,[a+b+c,a+b+c]+0.5,'color','k');
line([0,59]+0.5,[a+b+c+d,a+b+c+d]+0.5,'color','k');

box on;
set(gca,'XTick',[],'YTick',[])
box on;




%% Saccadic --> IBN (without function)
load ABDiPutativeSaccadic.mat
load ABD_ABDi_PutativeSaccadic.mat

SaccadicToIBN.IBNcellID = IBNall;

for i = 1:size(SaccadicAxons,2)
    SaccadicToIBN.SaccadicCellID(i) = SaccadicAxons(i);
    SaccadicToIBN.IBNcount(i,:) = histcounts(categorical(Sacc(i).Outputs),categorical(SaccadicToIBN.IBNcellID));
end


SaccadicToIBN.SacABDcounts = SaccadicToIBN.IBNcount(SaccadicProjectingToABDexclusively.Index,:);
SaccadicToIBN.SacABDicounts = SaccadicToIBN.IBNcount(SaccadicProjectingToABDiexclusively.Index,:);

SaccadicToIBN.AllABDcounts = sum(SaccadicToIBN.SacABDcounts,1);
SaccadicToIBN.AllABDicounts = sum(SaccadicToIBN.SacABDicounts,1);
 
 SaccadicToIBN.NormalizedCounts = (SaccadicToIBN.AllABDcounts -  SaccadicToIBN.AllABDicounts) ./ (SaccadicToIBN.AllABDcounts + SaccadicToIBN.AllABDicounts)

 figure;
 subplot(4,4,1)
 histogram(SaccadicToIBN.NormalizedCounts,-1:0.1:1,'FaceColor','k');
 box off;
 axis square;
 xlabel('(Sa-Sai)/(Sa+Sai)');
 
 SaccadicToIBN.Sac_ABDcellIDs = SaccadicToIBN.IBNcellID(SaccadicToIBN.NormalizedCounts>0);
 SaccadicToIBN.Sac_ABDicellIDs = SaccadicToIBN.IBNcellID(SaccadicToIBN.NormalizedCounts<0);
 
%  figure;
 %transform_swc_AV(SaccadicToIBN.Sac_ABDcellIDs,IbnABDcolor,[],true,true,'SaccIBNABD-Vel');
%  subplot(1,2,2)
%  transform_swc_AV(SaccadicToIBN.Sac_ABDicellIDs,IbnABDicolor,[],true,true,'SaccIBNABDi-Vel');



%% Connectivity

clear MatOrder
clear MatIndex
clear connMat


load AllCells.mat
load ConnMatrixPre.mat

load ABDVestPop.mat
load ABDiVestPop.mat
load ABDContraPop.mat
load ABDiContrapop.mat
load ABDRestPop.mat
load ABDiRestPop.mat
load ABDPutativeSaccadic.mat
load ABDiPutativeSaccadic.mat
load ABD_ABDi_PutativeSaccadic.mat

vestibularCellIds = [76631,76700,76656,76655,76660,76679,77393,77395,77375,...
     77815,77803,77807,80591,80579,77255,77697,77364,80223,81126,80760,81117,...
     80672,80622,80572,80569,81122,78426,81328,81575,81582,81870,82161,82165,81553,82169]; 
 
 MVNs = [76664 76784 76783 77394 77326 77772 77766 77756 77753 77767 77775 ...
     77780 80883 80247 80669 80698 80762 80748 80696 80761 80749 81070 81044 ...
     80991 80967 81008 80883 81045 80986 81023 78647 78619 80247 81141];
 
 lateralVSaccadic = [80163 80167 80177 80179 76688 80204 80206 80210 76682];
 
bushySaccadicLateral = [79067 79058 79069 79072 79055 79074 79062 79077 79085 ...
    79080 79086 79064 78576 78540 77146 78646 78542 78541 78543 80728]; 


%%
% sort Contra pop
tempA = [vertcat(ABDContrapop.cellIDs);vertcat(ABDiContrapop.cellIDs)];
tempA = unique(tempA);
ABDContraPop.rhombomereOrder = isRhombomere(tempA');
ABDContraPop.cellIDrhombomereOrdered = [tempA(logical(ABDContraPop.rhombomereOrder.r3));...
tempA(logical(ABDContraPop.rhombomereOrder.r4));tempA(logical(ABDContraPop.rhombomereOrder.r5));...
tempA(logical(ABDContraPop.rhombomereOrder.r6));tempA(logical(ABDContraPop.rhombomereOrder.r7))];
if size(ABDContraPop.cellIDrhombomereOrdered,1) ~= size(ABDContraPop.rhombomereOrder.cellID,1) 
    ABDContraPop.cellIDrhombomereOrdered(size(ABDContraPop.cellIDrhombomereOrdered,1)+1:size(ABDContraPop.rhombomereOrder.cellID,1)) = ...
        setdiff(tempA,ABDContraPop.cellIDrhombomereOrdered);
end
ABDContraPop.MotorDist = isMotor(tempA,df);
ABDContrPop.OSI = (sum(ABDContraPop.MotorDist(:,2:3),2)-sum(ABDContraPop.MotorDist(:,4:5),2)) ./ ...
    (sum(ABDContraPop.MotorDist(:,2:3),2)+sum(ABDContraPop.MotorDist(:,4:5),2));
[~,ABDContrPop.OSIsortOrder] = sort(ABDContrPop.OSI);
ABDContraPop.cellIDOSIrhombomereOrdered = tempA(flip(ABDContrPop.OSIsortOrder));
%ABDContraPop.cellIDOSIrhombomereOrdered = ABDContraPop.cellIDrhombomereOrdered(ABDContrPop.OSIsortOrder);

% sort Rest pop
tempB = [vertcat(ABDRestPop);vertcat(ABDiRestPop)];
clear ABDRestPop;
clear ABDiRestPop
tempB = unique(tempB);
ABDRestPop.cellIDs = tempB;
ABDRestPop.cellIDs = setdiff(ABDRestPop.cellIDs,ABD_ABDi_PutativeSaccadic.cellIDs);
ABDRestPop.MotorDist = isMotor(ABDRestPop.cellIDs,df);
ABDRestPop.OSI = (sum(ABDRestPop.MotorDist(:,2:3),2)-sum(ABDRestPop.MotorDist(:,4:5),2)) ./ ...
    (sum(ABDRestPop.MotorDist(:,2:3),2)+sum(ABDRestPop.MotorDist(:,4:5),2));
[~,ABDRestPop.OSIsortOrder] = sort(ABDRestPop.OSI);
ABDRestPop.cellIDOSIOrdered = ABDRestPop.cellIDs(flip(ABDRestPop.OSIsortOrder));


 %%
 
clear MatOrder
clear MatIndex
clear connMatCKT
clear connMatOrderedTotalPostSynapses;
clear ipsir78Int;
clear contrar78Int

ipsir78Int = [SaccadicToIntegrator.SaccABDIntCellIDs(isALX(SaccadicToIntegrator.SaccABDIntCellIDs));...
    SaccadicToIntegrator.SaccABDiIntCellIDs(isALX(SaccadicToIntegrator.SaccABDiIntCellIDs))];

contrar78Int = [SaccadicToIntegrator.SaccABDIntCellIDs(~isALX(SaccadicToIntegrator.SaccABDIntCellIDs));...
    SaccadicToIntegrator.SaccABDiIntCellIDs(~isALX(SaccadicToIntegrator.SaccABDiIntCellIDs))];

% PutativeSaccadic; Vestibular, IBNs, r456 Positon, Integrators, Motors, Contra
% Everythingelse
MatOrder = [ABDPutativeSaccadic.cellIDs;...
            ABDiPutativeSaccadic.cellIDs;...
            SaccadicToIBN.Sac_ABDcellIDs';...
            SaccadicToIBN.Sac_ABDicellIDs';...
            vestibularCellIds';...
            MVNs';...
            SaccadicProjectingToABDexclusively.cellID';...
            SaccadicProjectingToABDiexclusively.cellID';...
            SaccadicProjectingToABD_ABDi.cellID';...
            lateralVSaccadic';...
            ipsir78Int;...
            contrar78Int;...
            ABDr_CellIDs';ABDc_CellIDs';ABDIr_CellIDs';ABDIc_CellIDs';...
            %ABDContraPop.cellIDOSIrhombomereOrdered;...
            %ABDRestPop.cellIDOSIOrdered;...            
            ];
        
        
% Organize along the RC axis        
        
% MatOrderRC =   [ABDPutativeSaccadic.cellIDs;...
%             ABDiPutativeSaccadic.cellIDs;...
%             vestibularCellIds';...
%             MVNs';...
%             SaccadicToIBN.Sac_ABDcellIDs';...
%             SaccadicToIBN.Sac_ABDicellIDs';...
%             SaccadicProjectingToABDexclusively.cellID';...
%             SaccadicProjectingToABDiexclusively.cellID';...
%             SaccadicProjectingToABD_ABDi.cellID';...
%             lateralVSaccadic';...
%             ipsir78Int;...
%             contrar78Int;...
%             ABDr_CellIDs';ABDc_CellIDs';ABDIr_CellIDs';ABDIc_CellIDs';...
%             %ABDContraPop.cellIDOSIrhombomereOrdered;...
%             %ABDRestPop.cellIDOSIOrdered;...            
%             ];      
%         

 for i = 1:numel(MatOrder)
 temp = SynapticPartners(MatOrder(i),1,df);
connMatOrderedTotalPostSynapses(i,1) =  numel(temp);
connMatOrderedTotalReconPostSynapses(i,1) = numel(temp(temp<1e5));
clear temp;
end
        
for i = 1:numel(MatOrder)
    MatIndex(i) = find(AllCells == MatOrder(i),1);
end


endr23SaccABD = find(ABDPutativeSaccadic.cellIDs(end) == AllCells(MatIndex))+0.5;
endr23SaccABDi = find(ABDiPutativeSaccadic.cellIDs(end) == AllCells(MatIndex))+0.5;
endVestABD = find(vestibularCellIds(end) == AllCells(MatIndex))+0.5;
endMVN = find(MVNs(end) == AllCells(MatIndex))+0.5;
endIBNABD = find(SaccadicToIBN.Sac_ABDcellIDs(end) == AllCells(MatIndex))+0.5;
endIBNABDi = find(SaccadicToIBN.Sac_ABDicellIDs(end) == AllCells(MatIndex))+0.5;
endSacABD = find(SaccadicProjectingToABDexclusively.cellID(end) == AllCells(MatIndex))+0.5;
endSACABDi = find(SaccadicProjectingToABDiexclusively.cellID(end) == AllCells(MatIndex))+0.5;
endSACABD_ABDi = find(SaccadicProjectingToABD_ABDi.cellID(end) == AllCells(MatIndex))+0.5;
endSacVent = find(lateralVSaccadic(end) == AllCells(MatIndex))+0.5;
endSacIntABD = find(ipsir78Int(end) == AllCells(MatIndex))+0.5;
endSacIntABDi = find(contrar78Int(end) == AllCells(MatIndex))+0.5;

% endSacIntABD = find(SaccadicToIntegrator.SaccABDIntCellIDs(end) == AllCells(MatIndex))+0.5;
% endSacIntABDi = find(SaccadicToIntegrator.SaccABDiIntCellIDs(end) == AllCells(MatIndex))+0.5;

endABD = find(ABDc_CellIDs(end) == AllCells(MatIndex))+0.5;
endABDi = find(ABDIc_CellIDs(end) == AllCells(MatIndex))+0.5;

groupEnds= [endr23SaccABD,endr23SaccABDi,endVestABD,endMVN,endIBNABD,endIBNABDi,endSacABD,endSACABDi,endSACABD_ABDi,endSacVent,endSacIntABD,endSacIntABDi,endABD,endABDi];
groupIDs = {'r23SaccABD' 'r23SaccABDi' 'DON' 'MVN' 'IBN_ABD' 'IBN_ABDi' 'r456IntABD' 'r456IntABDi' 'r456IntABD_ABDi' 'r56VenInt' 'r78ALX' 'r78DBX' 'ABD' 'ABDi'};

%%
%subplot(2,2,1)
connMatCKT = ConnMatrixPre(MatIndex,MatIndex);
%cspy(connMat,'Colormap',colorcet('L2'),'Levels',255,'MarkerSize',10);
subplot(2,2,1);
%heatmap(connMat,'GridVisible','off','Colormap',colorcet('L4'));
cspy(connMatCKT,'Colormap',colorcet('R3','N',15),'Levels',15,'MarkerSize',7);
axis square;
% colormap(colorcet('L17','N',15,'reverse',1));
lighting phong;
material shiny;
set(gca, 'Xtick',[],'XTickLabel',[],'YTick',[],'YTickLabel',[]);
set(gca,'color',[0.90,0.9,0.90]);
colorbar;
box on;


hold on;
line([0,size(connMatCKT,1)],[endSacVent,endSacVent],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([0,size(connMatCKT,1)],[endr23SaccABD(1),endr23SaccABD(1)],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([0,size(connMatCKT,1)],[endr23SaccABDi,endr23SaccABDi],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([0,size(connMatCKT,1)],[endSacABD,endSacABD],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([0,size(connMatCKT,1)],[endSACABDi,endSACABDi],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([0,size(connMatCKT,1)],[endSACABD_ABDi,endSACABD_ABDi],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([0,size(connMatCKT,1)],[endVestABD(1),endVestABD(1)],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([0,size(connMatCKT,1)],[endMVN,endMVN],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([0,size(connMatCKT,1)],[endIBNABD,endIBNABD],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([0,size(connMatCKT,1)],[endIBNABDi,endIBNABDi],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([0,size(connMatCKT,1)],[endSacIntABD,endSacIntABD],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([0,size(connMatCKT,1)],[endSacIntABDi,endSacIntABDi],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([0,size(connMatCKT,1)],[endABD,endABD],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([0,size(connMatCKT,1)],[endABDi,endABDi],'color',[0.5,0.5,0.5],'lineWidth',0.5);


% 
line([endSacVent,endSacVent],[0,size(connMatCKT,1)],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([endr23SaccABD(1),endr23SaccABD(1)],[0,size(connMatCKT,1)],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([endr23SaccABDi,endr23SaccABDi],[0,size(connMatCKT,1)],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([endSacABD,endSacABD],[0,size(connMatCKT,1)],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([endSACABDi,endSACABDi],[0,size(connMatCKT,1)],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([endSACABD_ABDi,endSACABD_ABDi],[0,size(connMatCKT,1)],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([endVestABD(1),endVestABD(1)],[0,size(connMatCKT,1)],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([endMVN,endMVN],[0,size(connMatCKT,1)],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([endIBNABD,endIBNABD],[0,size(connMatCKT,1)],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([endIBNABDi,endIBNABDi],[0,size(connMatCKT,1)],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([endSacIntABD,endSacIntABD],[0,size(connMatCKT,1)],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([endSacIntABDi,endSacIntABDi],[0,size(connMatCKT,1)],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([endABD,endABD],[0,size(connMatCKT,1)],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([endABDi,endABDi],[0,size(connMatCKT,1)],'color',[0.5,0.5,0.5],'lineWidth',0.5);

subplot(2,2,2);
connMatNorm = bsxfun(@rdivide,connMatCKT,connMatOrderedTotalPostSynapses);
connMatNorm(connMatNorm>0.05) = 0.05;
cspy(connMatNorm ,'Colormap',colorcet('R3','N',15),'Levels',15,'MarkerSize',7);
caxis([0,0.05]);
axis square;
%colormap(colorcet('L17','N',15,'reverse',1));
lighting phong;
material shiny;
set(gca, 'Xtick',[],'XTickLabel',[],'YTick',[],'YTickLabel',[]);
set(gca,'color',[0.9,0.9,0.9]);
colorbar;
box on;


hold on;
line([0,size(connMatCKT,1)],[endSacVent,endSacVent],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([0,size(connMatCKT,1)],[endr23SaccABD(1),endr23SaccABD(1)],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([0,size(connMatCKT,1)],[endr23SaccABDi,endr23SaccABDi],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([0,size(connMatCKT,1)],[endSacABD,endSacABD],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([0,size(connMatCKT,1)],[endSACABDi,endSACABDi],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([0,size(connMatCKT,1)],[endSACABD_ABDi,endSACABD_ABDi],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([0,size(connMatCKT,1)],[endVestABD(1),endVestABD(1)],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([0,size(connMatCKT,1)],[endMVN,endMVN],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([0,size(connMatCKT,1)],[endIBNABD,endIBNABD],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([0,size(connMatCKT,1)],[endIBNABDi,endIBNABDi],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([0,size(connMatCKT,1)],[endSacIntABD,endSacIntABD],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([0,size(connMatCKT,1)],[endSacIntABDi,endSacIntABDi],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([0,size(connMatCKT,1)],[endABD,endABD],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([0,size(connMatCKT,1)],[endABDi,endABDi],'color',[0.5,0.5,0.5],'lineWidth',0.5);


% 
line([endSacVent,endSacVent],[0,size(connMatCKT,1)],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([endr23SaccABD(1),endr23SaccABD(1)],[0,size(connMatCKT,1)],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([endr23SaccABDi,endr23SaccABDi],[0,size(connMatCKT,1)],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([endSacABD,endSacABD],[0,size(connMatCKT,1)],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([endSACABDi,endSACABDi],[0,size(connMatCKT,1)],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([endSACABD_ABDi,endSACABD_ABDi],[0,size(connMatCKT,1)],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([endVestABD(1),endVestABD(1)],[0,size(connMatCKT,1)],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([endMVN,endMVN],[0,size(connMatCKT,1)],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([endIBNABD,endIBNABD],[0,size(connMatCKT,1)],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([endIBNABDi,endIBNABDi],[0,size(connMatCKT,1)],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([endSacIntABD,endSacIntABD],[0,size(connMatCKT,1)],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([endSacIntABDi,endSacIntABDi],[0,size(connMatCKT,1)],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([endABD,endABD],[0,size(connMatCKT,1)],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([endABDi,endABDi],[0,size(connMatCKT,1)],'color',[0.5,0.5,0.5],'lineWidth',0.5);

subplot(2,2,4)
barh(1:size(connMatNorm,1),sum(connMatNorm,2),'FaceColor','k');
set(gca, 'YDir','reverse','Xlim',[0,0.5],'Xtick',[0,0.25,0.5],'ylim',[0,size(connMatNorm,1)]);
box off;
daspect([1,160,1]);

%%
lineNumber = 313;
clear MatOrder
clear MatIndex
clear connMat

% PutativeSaccadic; Vestibular, IBNs, r456 Positon, Integrators, Motors, Contra
% Everythingelse
MatOrder = [ABDPutativeSaccadic.cellIDs;...
            ABDiPutativeSaccadic.cellIDs;...
            vestibularCellIds';...
            MVNs';...
            SaccadicToIBN.Sac_ABDcellIDs';...
            SaccadicToIBN.Sac_ABDicellIDs';...
            SaccadicProjectingToABDexclusively.cellID';...
            SaccadicProjectingToABDiexclusively.cellID';...
            SaccadicProjectingToABD_ABDi.cellID';...
            lateralVSaccadic';...
            SaccadicToIntegrator.SaccABDIntCellIDs;...
            SaccadicToIntegrator.SaccABDiIntCellIDs;...
            ABDr_CellIDs';ABDc_CellIDs';ABDIr_CellIDs';ABDIc_CellIDs';...
            ABDContraPop.cellIDOSIrhombomereOrdered;...
            ABDRestPop.cellIDOSIOrdered;...            
            ];

 for i = 1:numel(MatOrder)
 temp = SynapticPartners(MatOrder(i),1,df);
connMatOrderedTotalPostSynapses(i,1) =  numel(temp);
connMatOrderedTotalReconPostSynapses(i,1) = numel(temp(temp<1e5));
clear temp;
end
        
for i = 1:numel(MatOrder)
    MatIndex(i) = find(AllCells == MatOrder(i),1);
end

connMatCKT = ConnMatrixPre(MatIndex,MatIndex);
figure;
subplot(2,3,[1,2,4,5]);
cspy(connMatCKT,'Colormap',colorcet('L2','N',15,'reverse',1),'Levels',15,'MarkerSize',7,'Marker','s','filled');
axis square;
lighting phong;
material shiny;
set(gca, 'Xtick',[],'XTickLabel',[],'YTick',[],'YTickLabel',[]);
%set(gca,'color',[0.9,0.9,0.9]);
colorbar;
box on;

hold on;
line([0,size(connMatCKT,1)],[lineNumber,lineNumber],'color','k','lineWidth',0.5);
line([lineNumber,lineNumber],[0,size(connMatCKT,1)],'color','k','lineWidth',0.5);

subplot(2,3,[3,6]);
barh(1:size(connMatCKT,1),sum(connMatCKT,2),'FaceColor','k');
set(gca, 'YDir','reverse')%,'Xlim',[0,0.5],'Xtick',[0,0.25,0.5],'ylim',[0,size(connMatNorm,1)]);
box off;
%daspect([1,160,1]);



%% Vestibular projection patterns

Vestibular.DONcellIDs = vestibularCellIds;
Vestibular.MVNs = MVNs;

Vestibular.DONMotorCounts = isMotor(Vestibular.DONcellIDs',df);
Vestibular.MVNMotorCounts = isMotor(Vestibular.MVNs',df);

Vestibular.DONNormalized = ((Vestibular.DONMotorCounts(:,2)+Vestibular.DONMotorCounts(:,3)) - ((Vestibular.DONMotorCounts(:,4)+Vestibular.DONMotorCounts(:,5)))) ./ ...
    ((Vestibular.DONMotorCounts(:,2)+Vestibular.DONMotorCounts(:,3)) + (Vestibular.DONMotorCounts(:,4)+Vestibular.DONMotorCounts(:,5)));
    
    
Vestibular.MVNNormalized = ((Vestibular.MVNMotorCounts(:,2)+Vestibular.MVNMotorCounts(:,3)) - ((Vestibular.MVNMotorCounts(:,4)+Vestibular.MVNMotorCounts(:,5)))) ./ ...
    ((Vestibular.MVNMotorCounts(:,2)+Vestibular.MVNMotorCounts(:,3)) + (Vestibular.MVNMotorCounts(:,4)+Vestibular.MVNMotorCounts(:,5)));

figure;

subplot(4,4,1)
histogram(Vestibular.DONNormalized,-1:0.1:1,'FaceColor',lightMagenta);
hold on;
histogram(Vestibular.MVNNormalized,-1:0.1:1,'FaceColor',lightGreen);
axis square;
box off;

figure;
transform_swc_AV(Vestibular.DONcellIDs,colors(1,:),[],true,true);

%% ventral class of Sacc

lateralVSaccadic = [80163 80167 80177 80179 76688 80204 80206 80210 76682];
LateralVentralSaccadicMotorDist = isMotor(lateralVSaccadic',df);



bushySaccadicLateral = [79067 79058 79069 79072 79055 79074 79062 79077 79085 ...
    79080 79086 79064 78576 78540 77146 78646 78542 78546 78541 78543 80728]; 

bushySaccadicLateralMotorDist = isMotor(bushySaccadicLateral',df);

for i = 1:numel(bushySaccadicLateral)
    [a,~] = SynapticPartners(bushySaccadicLateral(i),1,df);
    bushySaccadicLateralSm(i) = sum(ismember(a,SaccadicProjectingToABDexclusively.cellID));
    bushySaccadicLateralSi(i) = sum(ismember(a,SaccadicProjectingToABDiexclusively.cellID));
    clear a;
end
