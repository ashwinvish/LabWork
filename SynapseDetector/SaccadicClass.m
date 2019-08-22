% Classes by anatomy
%clear;
addpath(genpath('/Users/ashwin/Documents/'));

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

%lateralVSaccadic = [80163 80167 80177 80179 76688 80204 80206 80210 76682];

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

confirmedDBX = [76199 76200 76182 76183 76185 76186 76189 76191 76188 ];
putativeDBX = [77582 77588 77595 77597 77598 77599 77605 79040 76399 76828 ...
    76829 78435 76826 76874 76289 76542 76832 76836 76838 76877 76883 76884 ...
    78400 78429 76383 76867 77594 76830 76876 78440 78438 78447 79045 78434 ...
    78450 76878 77683 77450 78448 ];

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

ABDvols = sortrows(ABDvols,2);
ABDvols(:,2) = ABDvols(:,2)./1e9;

for i = 1:size(ABDvols,1)
    ABDvols(i,3) = size(SynapticPartners(ABDvols(i,1),1,df),1);
    if (isExistReRoot(ABDvols(i,1)) == 1)
        tempTree = SwctoZbrian(ABDvols(i,1));
        ABDvols(i,4) = max(Pvec_tree(tempTree{1}));
        clear tempTree;
        ABDvols(i,5:7) = getOrigin(ABDvols(i,1));
    else
        ABDvols(i,4)= NaN;  
    end
    %ABDvols(i,5:8) = getOrigin(ABDvols(i,1));
end

ABDr_vols = ABDvols(find(ismember(ABDvols(:,1),ABDr_CellIDs)),2);
ABDr_TotalSynapses = ABDvols(find(ismember(ABDvols(:,1),ABDr_CellIDs)),3);

ABDc_vols = ABDvols(find(ismember(ABDvols(:,1),ABDc_CellIDs)),2);
ABDc_TotalSynapses = ABDvols(find(ismember(ABDvols(:,1),ABDc_CellIDs)),3);

ABDIr_vols = ABDvols(find(ismember(ABDvols(:,1),ABDIr_CellIDs)),2);
ABDIr_TotalSynapses = ABDvols(find(ismember(ABDvols(:,1),ABDIr_CellIDs)),3);

ABDIc_vols = ABDvols(find(ismember(ABDvols(:,1),ABDIc_CellIDs)),2);
ABDIc_TotalSynapses = ABDvols(find(ismember(ABDvols(:,1),ABDIc_CellIDs)),3);

figure;
subplot(4,4,1)
scatter(ABDvols(:,2),ABDvols(:,3),'ko');
xlabel('ABD vol (um^3)');
ylabel('Total synapses');
axis square
subplot(4,4,3)
scatter(ABDvols(:,2),ABDvols(:,4),'ko');
xlabel('ABD vol (um^3)');
ylabel('ABD path length (um)');
axis square
subplot(4,4,6)
scatter(ABDvols(:,3),ABDvols(:,4),'ko');
scatter(ABDvols(:,4),ABDvols(:,3),'ko');
xlabel('path length (um)');
ylabel('Total synapses');
axis square
subplot(4,4,8)
scatter3(ABDvols(:,2),ABDvols(:,4),ABDvols(:,3),'ko');
xlabel('ABD vol (um^3)');
ylabel('Pathlength (um)');
zlabel('Total Synapses');
axis square
box on

%%
 for i = 1:numel(SaccadicAxons)
     if ~isExistReRoot(SaccadicAxons(i))==0
     Sacc(i) = InputsByClass(SaccadicAxons(i),df,2);
     end  
 end
%%

SaccadicMotorDistribution.allCounts = vertcat(Sacc.MotorDist);
noSynapsesOnMotor = find(sum(SaccadicMotorDistribution.allCounts(:,2:5),2)==0);
SaccadicMotorDistribution.allCounts(noSynapsesOnMotor,:) = NaN;
SaccadicMotorDistribution.ABD = sum(SaccadicMotorDistribution.allCounts(:,2:3),2);
SaccadicMotorDistribution.ABDi = sum(SaccadicMotorDistribution.allCounts(:,4:5),2);
SaccadicMotorDistribution.normalizedMotorCounts(:,1) = SaccadicMotorDistribution.allCounts(:,1);
SaccadicMotorDistribution.normalizedMotorCounts(:,2) = (SaccadicMotorDistribution.ABD-SaccadicMotorDistribution.ABDi)./ (SaccadicMotorDistribution.ABD+SaccadicMotorDistribution.ABDi);
SaccadicMotorDistribution.normalized(:,2) = SaccadicMotorDistribution.normalizedMotorCounts(~isnan(SaccadicMotorDistribution.normalizedMotorCounts(:,2)),2);
SaccadicMotorDistribution.normalized(:,1) = SaccadicMotorDistribution.normalizedMotorCounts(~isnan(SaccadicMotorDistribution.normalizedMotorCounts(:,2)),1); 



[~,sortSac] = sortrows(SaccadicMotorDistribution.normalized,2);
cmap = colorcet('D2','N',numel(unique(SaccadicMotorDistribution.normalized(:,2))),'reverse',1);
cmapNew = [repmat(cmap(1,:),31,1);cmap;repmat(cmap(end,:),19,1)];

transform_swc_AV(SaccadicMotorDistribution.normalized(sortSac,1),cmapNew,[],true,false);
colormap(gca,cmapNew);
colorbar(gca);

SaccadicProjectingToABDexclusively = SaccadicMotorDistribution.normalizedMotorCounts(find(SaccadicMotorDistribution.normalizedMotorCounts(:,2) == 1));
SaccadicProjectingToABDiexclusively = SaccadicMotorDistribution.normalizedMotorCounts(find(SaccadicMotorDistribution.normalizedMotorCounts(:,2) ==-1));
SaccadicProjectingToABD_ABDi = SaccadicMotorDistribution.normalizedMotorCounts(find(SaccadicMotorDistribution.normalizedMotorCounts(:,2)>-0.95 & SaccadicMotorDistribution.normalizedMotorCounts(:,2)<0.95 ))
 
figure;
subplot(4,4,1)
scatter(SaccadicMotorDistribution.ABD,SaccadicMotorDistribution.ABDi,'ko','filled');
axis square;
set(gca,'XLim',[0,150],'ylim',[0,150]);
xlabel('Synapses on M');
ylabel('Synapses on I');
offsetAxes(gca);

subplot(4,4,2)
clear temp;
temp = isMotor(SaccadicProjectingToABDexclusively,df);
tempSum = sum(temp(:,2:3),2);
histogram(tempSum,0:20:150,'DisplayStyle','stairs','EdgeColor',SaccABDcolor);
hold on;
clear temp;
temp = isMotor(SaccadicProjectingToABDiexclusively,df);
tempSum = sum(temp(:,2:3),2);
histogram(tempSum,0:20:150,'DisplayStyle','stairs','EdgeColor',SaccABDicolor);
clear temp;
temp = isMotor(SaccadicProjectingToABD_ABDi,df);
tempSum = sum(temp(:,2:3),2);
histogram(tempSum,0:20:150,'DisplayStyle','stairs','EdgeColor',[0.8,0.8,0.8]);
box off;
axis square;
legend('Sm','Si','both');
offsetAxes(gca);



ix1 = 1;
ix2 = 1;
for i = 1:size(SaccadicMotorDistribution.allCounts,1)
    if SaccadicMotorDistribution.normalizedMotorCounts(i,2)>0.5
    ABDSaccadicpop(ix1).cellIDs = SaccadicMotorDistribution.allCounts(i,1);
    ABDSaccadicpop(ix1).Origin = Sacc(find(ABDSaccadicpop(ix1).cellIDs==SaccadicAxons)).Origin;
    ABDSaccadicpop(ix1).MotorCounts = SaccadicMotorDistribution.normalizedMotorCounts(i,2);
    ix1 = ix1+1;
    elseif SaccadicMotorDistribution.normalizedMotorCounts(i,2)<-0.5
     ABDiSaccadicpop(ix2).cellIDs = SaccadicMotorDistribution.allCounts(i,1);
     ABDiSaccadicpop(ix2).Origin = Sacc(find(ABDiSaccadicpop(ix2).cellIDs==SaccadicAxons)).Origin;
     ABDiSaccadicpop(ix2).MotorCounts = SaccadicMotorDistribution.normalizedMotorCounts(i,2);
     ix2 = ix2+1;
    end
end

 SaccadicProjectingToABD = [ABDSaccadicpop.cellIDs];
 SaccadicProjectingToABDi = [ABDiSaccadicpop.cellIDs];
 
 
 % make individual plots
 
 for i = 1:size(SaccadicProjectingToABD_ABDi,1)
     %subplot(4,4,i)
     figure('units','normalized','outerposition',[0 0 1 1]);
     transform_swc_AV(SaccadicProjectingToABD_ABDi(i),[0,0,0],[],true,false);
     fname = sprintf('/Users/ashwin/Desktop/%5s.png',SaccadicProjectingToABD_ABDi(i))
     export_fig(fname);
     %hold on;
     %transform_swc_AV([ABDc_CellIDs,ABDr_CellIDs],ABDcolor,[],false,false);
     %transform_swc_AV([ABDIc_CellIDs,ABDIr_CellIDs],ABDicolor,[],false,false);
 end
 

%%
figure;
subplot(1,2,1)
[~,SaccadicABDindex]= sort([ABDSaccadicpop.MotorCounts]); % 0-->1
SaccadicABDindex = flip(SaccadicABDindex); % 1-->0
cmapReds = colorcet('D1','N', 2*length(SaccadicABDindex));
%transform_swc_AV([ABDSaccadicpop(SaccadicABDindex).cellIDs],cmapReds,[],true,false);
transform_swc_AV(SaccadicProjectingToABDexclusively,SaccABDcolor,[],true,false);

colormap(gca,cmapReds(length(SaccadicABDindex):end,:))
colorbar(gca);
subplot(1,2,2)
clear index;
[~,SaccadicABDiindex]= sort([ABDiSaccadicpop.MotorCounts]); % -1--0
SaccadicABDiindex = flip(SaccadicABDiindex); % 0-->-1
%SaccadicABDiindex = flip(SaccadicABDiindex,2);
%cmapBlues = colorcet('D1','N', 2*length(SaccadicABDiindex));
%transform_swc_AV([ABDiSaccadicpop(SaccadicABDiindex).cellIDs],colors(2,:),[],true,false);
transform_swc_AV(SaccadicProjectingToABDiexclusively,SaccABDicolor,[],true,false);

%colormap(gca,flipud(cmapBlues(1:length(SaccadicABDiindex),:)));
%colorbar(gca);



H = hdf5read('/Users/ashwin/Documents/LabWork/SynapseDetector/contact_area.h5','/matrix');
HLables = hdf5read('/Users/ashwin/Documents/LabWork/SynapseDetector/contact_area.h5','/neuron_id_list');
Horder = ismember(HLables,SaccadicAxons);
Hsub = H(Horder(1:find(HLables==80315)),:); % till location of last saccadic axon
SaccadicABDcontactArea = sum(Hsub(:,find(HLables == 82140):find(HLables == 77296)),2);
SaccadicABDicontactArea = sum(Hsub(:,find(HLables == 78553):end),2);

NormalizedSaccadicContactArea = (SaccadicABDcontactArea-SaccadicABDicontactArea)./...
                                (SaccadicABDcontactArea+SaccadicABDicontactArea);
figure;
subplot(4,4,1)
histogram(SaccadicMotorDistribution.normalizedMotorCounts(:,2),-1:0.1:1,'FaceColor','k','FaceAlpha',0.5);
hold on;
line([0.5,0.5],[0,40],'color','k','LineStyle',':');
line([-0.5,-0.5],[0,40],'color','k','LineStyle',':');
%histogram(NormalizedSaccadicContactArea,-1:0.1:1,'FaceColor','r','FaceAlpha',0.5);
xlabel('(A-Ai)/A+Ai)');
axis square;
box off;

figure; 
subplot(4,4,2)
scatter(SaccadicABDcontactArea,SaccadicABDicontactArea,25,'ko','filled');
                            
location = [vertcat(ABDSaccadicpop.Origin);vertcat(ABDiSaccadicpop.Origin)];
groups = [ones(size(vertcat(ABDSaccadicpop.Origin),1),1);2*ones(size(vertcat(ABDiSaccadicpop.Origin),1),1)];

figure;
subplot(4,4,[3,4])
scatterhist(location(:,1),location(:,2),'Group',groups,'color',[lightRed;lightBlue]);
set(gca,'YDir','reverse');


 %% potential Synapses analysis
AllABDPreSynapticLocaions = [vertcat(ABDr.PreSynCoordsTransformed); vertcat(ABDc.PreSynCoordsTransformed)];
AllABDiPreSynapticLoactions = [vertcat(ABDIr.PreSynCoordsTransformed); vertcat(ABDIc.PreSynCoordsTransformed)];
Saccadic.ABDpotentialSynapses = zeros(numel(SaccadicAxons),2);
Saccadic.ABDipotentialSynapses = zeros(numel(SaccadicAxons),2);


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
        
    end
end
 
SaccadeABDexclusiveIndex = find(ismember(SaccadicAxons,SaccadicProjectingToABDexclusively));
SaccadeABDiexclusiveIndex = find(ismember(SaccadicAxons,SaccadicProjectingToABDiexclusively));
Saccade_ABD_ABDi_Index = find(ismember(SaccadicAxons, SaccadicProjectingToABD_ABDi));

ABDtruePotentialRatio = Saccadic.PotentialABDSynapses(:,2)./Saccadic.PotentialABDSynapses(:,1) ;
ABDiturePoentialRatio = Saccadic.PotentialABDiSynapses(:,2)./ Saccadic.PotentialABDiSynapses(:,1);

figure;
subplot(4,4,1)
scatter(Saccadic.PotentialABDSynapses(:,2),Saccadic.PotentialABDSynapses(:,1),20,SaccABDcolor,'filled');
hold on;
scatter(Saccadic.PotentialABDiSynapses(:,2),Saccadic.PotentialABDiSynapses(:,1),20,SaccABDicolor,'filled',...
    'Marker','s');
daspect([1,1,1]);
xlabel('Number of synapses');
ylabel('Potential Synapses');


subplot(4,4,2)
scatter(ABDtruePotentialRatio,ABDiturePoentialRatio,30,'ko','filled');
line([0,1],[0,1],'color','k','LineStyle','-');
line([0.5,0.5],[0,1],'color',[0.8,0.8,0.8],'LineStyle',':');
line([0,1],[0.5,0.5],'color',[0.8,0.8,0.8],'LineStyle',':');
axis square;
box off;
xlabel('synapses/potential onto M');
ylabel('synapses/potential onto I');
offsetAxes(gca);


subplot(4,4,3)
jitterValuesY = 2*(rand(size(SaccadeABDexclusiveIndex))-0.05)*0.05;   % +/-jitterAmount ma
 scatter(ABDtruePotentialRatio(SaccadeABDexclusiveIndex),ABDiturePoentialRatio(SaccadeABDexclusiveIndex)+jitterValuesY',...
     30,SaccABDcolor,'filled','MarkerEdgeColor',[0.5,0.5,0.5]);
 hold on;
 scatter(ABDtruePotentialRatio(SaccadeABDiexclusiveIndex), ABDiturePoentialRatio(SaccadeABDiexclusiveIndex),...
    30,SaccABDicolor,'filled','MarkerEdgeColor',[0.5,0.5,0.5],'jitter','on', 'jitterAmount',0.05); 
line([0,1],[0,1],'color','k','LineStyle','-');
line([0.5,0.5],[0,1],'color',[0.8,0.8,0.8],'LineStyle',':');
line([0,1],[0.5,0.5],'color',[0.8,0.8,0.8],'LineStyle',':');
axis square;
box off;
xlabel('synapses/potential onto M');
ylabel('synapses/potential onto I');
offsetAxes(gca);


subplot(4,4,5)
scatter(Saccadic.PotentialABDSynapses(SaccadeABDexclusiveIndex,2), Saccadic.PotentialABDiSynapses(SaccadeABDexclusiveIndex,1),30,SaccABDcolor,'filled');
hold on;
scatter(Saccadic.PotentialABDSynapses(Saccade_ABD_ABDi_Index,2), Saccadic.PotentialABDiSynapses(Saccade_ABD_ABDi_Index,1),30,'k','filled');
axis square;
box off;
xlabel('Actual Synapses on M');
ylabel('potential Synapses to I');
offsetAxes(gca);

subplot(4,4,6)
scatter(Saccadic.PotentialABDiSynapses(SaccadeABDiexclusiveIndex,2), Saccadic.PotentialABDSynapses(SaccadeABDiexclusiveIndex,1),30,SaccABDicolor,'filled');
hold on;
scatter(Saccadic.PotentialABDiSynapses(Saccade_ABD_ABDi_Index,2), Saccadic.PotentialABDSynapses(Saccade_ABD_ABDi_Index,1),30,'k','filled');
axis square;
box off;
xlabel('Actual Synapses on I');
ylabel('potential Synapses to M');
offsetAxes(gca);


%% Saccadic Motor partners

for i = 1:numel([ABDSaccadicpop.cellIDs])
  ABDSaccadicpop(i).MotorNeuronCounts = isPostSynapseMotor(ABDSaccadicpop(i).cellIDs,df); 
end

for i = 1:numel([ABDiSaccadicpop.cellIDs])
  ABDiSaccadicpop(i).MotorNeuronCounts = isPostSynapseMotor(ABDiSaccadicpop(i).cellIDs,df); 
end

for i = 1:numel([ABDSaccadicpop.cellIDs])
    ABDrSaccadicpopMotorConn(i,:) = ABDSaccadicpop(SaccadicABDindex(i)).MotorNeuronCounts(1,1:14)./ABDr_TotalSynapses';
    ABDcSaccadicpopMotorConn(i,:) = ABDSaccadicpop(SaccadicABDindex(i)).MotorNeuronCounts(2,1:17)./ABDc_TotalSynapses';
    ABDIrSaccadicpopMotorConn(i,:) = ABDSaccadicpop(SaccadicABDindex(i)).MotorNeuronCounts(3,1:11)./ABDIr_TotalSynapses';
    ABDIcSaccadicpopMotorConn(i,:) = ABDSaccadicpop(SaccadicABDindex(i)).MotorNeuronCounts(4,1:10)./ABDIc_TotalSynapses';
end

for i = 1:numel([ABDiSaccadicpop.cellIDs])
    ABDrSaccadicpopMotorConn(i+numel([ABDSaccadicpop.cellIDs]),:) = ABDiSaccadicpop(SaccadicABDiindex(i)).MotorNeuronCounts(1,1:14)./ABDr_TotalSynapses';
    ABDcSaccadicpopMotorConn(i+numel([ABDSaccadicpop.cellIDs]),:) = ABDiSaccadicpop(SaccadicABDiindex(i)).MotorNeuronCounts(2,1:17)./ABDc_TotalSynapses';
    ABDIrSaccadicpopMotorConn(i+numel([ABDSaccadicpop.cellIDs]),:) = ABDiSaccadicpop(SaccadicABDiindex(i)).MotorNeuronCounts(3,1:11)./ABDIr_TotalSynapses';
    ABDIcSaccadicpopMotorConn(i+numel([ABDSaccadicpop.cellIDs]),:) = ABDiSaccadicpop(SaccadicABDiindex(i)).MotorNeuronCounts(4,1:10)./ABDIc_TotalSynapses';
end


subplot(4,4,[1,5,9])
heatmap(ABDrSaccadicpopMotorConn);

subplot(4,4,[2,6,10])
heatmap(ABDcSaccadicpopMotorConn);


subplot(4,4,[3,7,11])
heatmap(ABDIrSaccadicpopMotorConn);


subplot(4,4,[4,8,12])
heatmap(ABDIcSaccadicpopMotorConn);


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

for i = 1:numel(ABDr_CellIDs)
    temp1 = ismember(ABDr(i).Inputs,vertcat([ABDSaccadicpop.cellIDs]',[ABDiSaccadicpop.cellIDs]'));
    ABDSaccPathLengthABDr{i} = ABDr(i).PathLength(temp1)'./max(Pvec_tree(ABDr(i).Tree{1}));
    clear temp1;
    ABDr_origin(i,:) = getOrigin(ABDr_CellIDs(i));
end


for i = 1:numel(ABDc_CellIDs)
    if ~isempty(ABDc(i).Tree)
    temp1 = ismember(ABDc(i).Inputs,vertcat([ABDSaccadicpop.cellIDs]',[ABDiSaccadicpop.cellIDs]'));
    ABDSaccPathLengthABDc{i} = ABDc(i).PathLength(temp1)'./max(Pvec_tree(ABDc(i).Tree{1}));
    clear temp1;
    end
    ABDc_origin(i,:) = getOrigin(ABDc_CellIDs(i));
end

for i = 1:numel(ABDIr_CellIDs)
    temp1 = ismember(ABDIr(i).Inputs,vertcat([ABDiSaccadicpop.cellIDs]',[ABDSaccadicpop.cellIDs]'));
    ABDSaccPathLengthABDIr{i} = ABDIr(i).PathLength(temp1)'./max(Pvec_tree(ABDIr(i).Tree{1}));
    clear temp1;
    ABDIr_origin(i,:) = getOrigin(ABDIr_CellIDs(i));
end


for i = 1:numel(ABDIc_CellIDs)
    temp1 = ismember(ABDIc(i).Inputs,vertcat([ABDiSaccadicpop.cellIDs]',[ABDSaccadicpop.cellIDs]'));
    ABDSaccPathLengthABDIc{i} = ABDIc(i).PathLength(temp1)'./max(Pvec_tree(ABDIc(i).Tree{1}));
    clear temp1;
    ABCIc_origin(i,:) = getOrigin(ABDIc_CellIDs(i));
end


for i = 1:numel(ABDr_CellIDs)
   ABDrSaccadicGradient(i,:) =  histcounts(ABDSaccPathLengthABDr{i},0:0.1:1);
   ABDrSaccadicGradient(i,:) = ABDrSaccadicGradient(i,:)/max(ABDrSaccadicGradient(i,:));
end

for i = 1:numel(ABDc_CellIDs)
   ABDcSaccadicGradient(i,:) =  histcounts(ABDSaccPathLengthABDc{i},0:0.1:1);
    ABDcSaccadicGradient(i,:) =  ABDcSaccadicGradient(i,:)/max( ABDcSaccadicGradient(i,:));
end

for i = 1:numel(ABDIr_CellIDs)
   ABDIrSaccadicGradient(i,:) =  histcounts(ABDSaccPathLengthABDIr{i},0:0.1:1);
   ABDIrSaccadicGradient(i,:) = ABDIrSaccadicGradient(i,:)/max(ABDIrSaccadicGradient(i,:));
end

for i = 1:numel(ABDIc_CellIDs)
   ABDIcSaccadicGradient(i,:) =  histcounts(ABDSaccPathLengthABDIc{i},0:0.1:1);
    ABDIcSaccadicGradient(i,:) =  ABDIcSaccadicGradient(i,:)/max(ABDIcSaccadicGradient(i,:));
end

figure;
subplot(4,4,[1,5,9,13])
imagesc(ABDrSaccadicGradient);
set(gca,'XTickLabel',[0.2,0.4,0.6,0.8,1]);
colormap(gca,colorcet('L8'));
colorbar(gca);

subplot(4,4,2)
histogram(ABDSaccPathLengthABDr{2},0:0.1:1,'Normalization','probability');
title('row:2')
box off;
subplot(4,4,6)
shadedErrorBar(0:0.1:0.9,mean(ABDrSaccadicGradient(1:7,:)),std(ABDrSaccadicGradient(1:7,:))/sqrt(7),'lineprops',{'Color','k'});
hold on;
shadedErrorBar(0:0.1:0.9,mean(ABDrSaccadicGradient(8:14,:)),std(ABDrSaccadicGradient(8:14,:))/sqrt(7),'lineprops',{'Color','r'});
box off;
subplot(4,4,14)
histogram(ABDSaccPathLengthABDr{13},0:0.1:1,'Normalization','probability');
title('row:13');
box off;

subplot(4,4,[3,7,11,15])
imagesc(ABDcSaccadicGradient);
set(gca,'XTickLabel',[0.2,0.4,0.6,0.8,1]);
colormap(gca,colorcet('L8'));
colorbar(gca);

subplot(4,4,4)
histogram(ABDSaccPathLengthABDc{2},0:0.1:1,'Normalization','probability');
title('row:2');
box off;
subplot(4,4,8)
shadedErrorBar(0:0.1:0.9,nanmean(ABDcSaccadicGradient(1:9,:)),nanstd(ABDcSaccadicGradient(1:9,:))/sqrt(8),'lineprops',{'Color','k'});
hold on;
shadedErrorBar(0:0.1:0.9,nanmean(ABDcSaccadicGradient(10:17,:)),nanstd(ABDcSaccadicGradient(10:17,:))/sqrt(7),'lineprops',{'Color','r'});
box off;
subplot(4,4,16)
histogram(ABDSaccPathLengthABDc{15},0:0.1:1,'Normalization','probability');
title('row:15');
box off;

figure;

subplot(4,4,[1,5,9,13])
imagesc(ABDIrSaccadicGradient);
set(gca,'XTickLabel',[0.2,0.4,0.6,0.8,1]);
colormap(gca,colorcet('L8'));
colorbar(gca);
subplot(4,4,2)
histogram(ABDSaccPathLengthABDIr{3},0:0.1:1,'Normalization','probability');
title('row:3');
box off;
subplot(4,4,6)
shadedErrorBar(0:0.1:0.9,mean(ABDIrSaccadicGradient(1:5,:)),std(ABDIrSaccadicGradient(1:5,:))/sqrt(5),'lineprops',{'Color','k'});
hold on;
shadedErrorBar(0:0.1:0.9,mean(ABDIrSaccadicGradient(6:11,:)),std(ABDIrSaccadicGradient(6:11,:))/sqrt(6),'lineprops',{'Color','r'});
box off;
subplot(4,4,14)
histogram(ABDSaccPathLengthABDIr{10},0:0.1:1,'Normalization','probability');
title('row:10');
box off;

subplot(4,4,[3,7,11,15])
imagesc(ABDIcSaccadicGradient);
set(gca,'XTickLabel',[0.2,0.4,0.6,0.8,1]);
colormap(gca,colorcet('L8'));
colorbar(gca);

subplot(4,4,4)
histogram(ABDSaccPathLengthABDIc{3},0:0.1:1,'Normalization','probability');
title('row:3');
box off;
subplot(4,4,8)
shadedErrorBar(0:0.1:0.9,nanmean(ABDIcSaccadicGradient(1:5,:)),nanstd(ABDIcSaccadicGradient(1:5,:))/sqrt(5),'lineprops',{'Color','k'});
hold on;
shadedErrorBar(0:0.1:0.9,nanmean(ABDIcSaccadicGradient(6:10,:)),nanstd(ABDIcSaccadicGradient(6:10,:))/sqrt(4),'lineprops',{'Color','r'});
box off;
subplot(4,4,16)
histogram(ABDSaccPathLengthABDIc{10},0:0.1:1,'Normalization','probability');
title('row:10');
box off;
%% plot integtators that are post synaptic to saccadic classes

SaccadicToIntegrator.IntegratorCellIDs = unique(vertcat(Sacc.Integrator));
SaccadicToABD.IntegratorIDs = [];
for i = 1:numel(SaccadicProjectingToABD)
    t = find([Sacc.cellID] == SaccadicProjectingToABD(i));
    SaccadicToABD.IntegratorIDs = [SaccadicToABD.IntegratorIDs;Sacc(t).Integrator];
    clear t;
end

SaccadicToABDi.IntegratorIDs = [];
for i = 1:numel(SaccadicProjectingToABDi)
    t = find([Sacc.cellID] == SaccadicProjectingToABDi(i));
    SaccadicToABDi.IntegratorIDs = [SaccadicToABDi.IntegratorIDs;Sacc(t).Integrator];
    clear t;
end

figure;
subplot(1,2,1) 
transform_swc_AV(unique(SaccadicToABD.IntegratorIDs),leadColor,[],true,false);
%colormap(gca,cmapReds(length(SaccadicABDindex):end,:))
%colorbar(gca);

subplot(1,2,2)
transform_swc_AV(unique(SaccadicToABDi.IntegratorIDs),lagColor,[],true,false);
%colormap(gca,flipud(cmapBlues(1:length(SaccadicABDiindex),:)));
%colorbar(gca);



%% Saccadic --> Integrator (without function)

SaccadicToIntegrator.IntegratorCellIDs = unique(vertcat(Sacc.Integrator));
   
for i = 1:size(SaccadicAxons,2)
    SaccadicToIntegrator.SaccadicCellID(i) = SaccadicAxons(i);
    SaccadicToIntegrator.IntegratorCount(i,:) = histcounts(categorical(Sacc(i).Integrator),categorical(SaccadicToIntegrator.IntegratorCellIDs));
end

SaccadicToIntegrator.MotorDistribution = isMotor(SaccadicToIntegrator.SaccadicCellID',df);

% find the saccadic neurons that project to ABD and ABDi
[~,SaccadicToIntegrator.ABDorder] = intersect(SaccadicToIntegrator.SaccadicCellID,SaccadicProjectingToABD,'stable');
[~,SaccadicToIntegrator.ABDiorder] = intersect(SaccadicToIntegrator.SaccadicCellID,SaccadicProjectingToABDi,'stable');

% ID of Saccadic axons that project to ABD and ABDi
SaccadicToIntegrator.ABDSaccadicCellIDs = SaccadicToIntegrator.SaccadicCellID(SaccadicToIntegrator.ABDorder);
%SaccadicToIntegrator.ABDSaccadicNormCount = SaccadicMotorDistribution.normalizedMotorCounts(SaccadicToIntegrator.ABDorder);
SaccadicToIntegrator.ABDiSaccadicCellIDs = SaccadicToIntegrator.SaccadicCellID(SaccadicToIntegrator.ABDiorder);
%SaccadicToIntegrator.ABDiSaccadicNormCount = SaccadicMotorDistribution.normalizedMotorCounts(SaccadicToIntegrator.ABDiorder);

% get counts onto integrators from axon projecting to ABD and ABDi
SaccadicToIntegrator.ABDIntCount = sum(SaccadicToIntegrator.IntegratorCount(SaccadicToIntegrator.ABDorder,:),1);
SaccadicToIntegrator.ABDiIntCount = sum(SaccadicToIntegrator.IntegratorCount(SaccadicToIntegrator.ABDiorder,:),1);
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

subplot(4,4,2)
histogram(SaccadicToIntegrator.ABDIntNormalizedCount(ismember(SaccadicToIntegrator.IntegratorCellIDs,allALX)),-1:0.1:1,'FaceColor','k')
box off;
xlabel('(Sa-Sai)/(Sa+Sai)');
axis square;

subplot(4,4,3)
histogram(SaccadicToIntegrator.ABDIntNormalizedCount(ismember(SaccadicToIntegrator.IntegratorCellIDs,allDBX)),-1:0.1:1,'FaceColor','k')
box off;
xlabel('(Sa-Sai)/(Sa+Sai)');
axis square;

figure;
subplot(1,2,1)
SaccadicToIntegrator.SaccABDIntOrder = find(SaccadicToIntegrator.ABDIntNormalizedCount>0);
SaccadicToIntegrator.SaccABDIntCellIDs = SaccadicToIntegrator.IntegratorCellIDs(SaccadicToIntegrator.SaccABDIntOrder);
%cmapReds = colorcet('D1','N', 2*length(SaccadicToIntegrator.SaccABDIntOrder));
%plotOrder = 
transform_swc_AV(SaccadicToIntegrator.SaccABDIntCellIDs,leadColor,[],true,false);
%colormap(gca,cmapReds(length(SaccadicABDindex):end,:))
%colorbar(gca);

subplot(1,2,2)
SaccadicToIntegrator.SaccABDiIntOrder = find(SaccadicToIntegrator.ABDIntNormalizedCount<0);
SaccadicToIntegrator.SaccABDiIntCellIDs = SaccadicToIntegrator.IntegratorCellIDs(SaccadicToIntegrator.SaccABDiIntOrder);
cmapBlues = colorcet('D1','N', 2*length(SaccadicABDiindex));
transform_swc_AV(SaccadicToIntegrator.SaccABDiIntCellIDs,lagColor,[],true,false);
%colormap(gca,flipud(cmapBlues(1:length(SaccadicABDiindex),:)));
%colorbar(gca);

figure;




figure;
imagesc(SaccadicToIntegrator.IntegratorCount([SaccadicToIntegrator.ABDorder;SaccadicToIntegrator.ABDiorder],[SaccadicToIntegrator.SaccABDIntOrder';SaccadicToIntegrator.SaccABDiIntOrder']))
colormap(colorcet('L8'));
line([0,size(SaccadicToIntegrator.IntegratorCellIDs,1)],[size(SaccadicToIntegrator.ABDorder,1),size(SaccadicToIntegrator.ABDorder,1)],'color','k');
line([size(SaccadicToIntegrator.SaccABDIntOrder,2),size(SaccadicToIntegrator.SaccABDIntOrder,2)],[0,size(SaccadicToIntegrator.SaccadicCellID,2)],'color','k');

%% Do integrators obey ABD order or Sac-->ABD order?

allIntegrators = [confirmedALX,putativeALX,confirmedDBX,putativeDBX,confirmedBARHL,putativeBARHL];
Integrators.motorDistribution = isMotor(allIntegrators',df);
Integrators.NormalizedCount = (sum(Integrators.motorDistribution(:,2:3),2)-sum(Integrators.motorDistribution(:,4:5),2))./...
                               (sum(Integrators.motorDistribution(:,2:3),2)+sum(Integrators.motorDistribution(:,4:5),2)) ;

[~,allIntegratorOrder] = ismember(allIntegrators,SaccadicToIntegrator.IntegratorCellIDs);
allIntegratorOrder(allIntegratorOrder==0) = [];

ALXids = 68;
DBXids = 68+48;
BarhlIds = 68+48+9;

IntegratorColors = [repmat(colors(1,:),ALXids,1);repmat(colors(2,:),DBXids,1);repmat(colors(3,:),BarhlIds,1)];

figure;
subplot(4,4,1);
scatter(Integrators.NormalizedCount(allIntegratorOrder),SaccadicToIntegrator.ABDIntNormalizedCount,20,IntegratorColors(allIntegratorOrder',:),'filled');
offsetAxes(gca);
xlabel('$$\frac{(Ia-Ii)}{(Ia+Ii)}$$','Interpreter','latex');
ylabel('Sa-Si)/(Sa+Si)');
axis square;
box off;

%% Saccadic --> IBN (without function)

SaccadicToIBN.IBNcellID = IBNall;

for i = 1:size(SaccadicAxons,2)
    SaccadicToIBN.SaccadicCellID(i) = SaccadicAxons(i);
    SaccadicToIBN.IBNcount(i,:) = histcounts(categorical(Sacc(i).Outputs),categorical(SaccadicToIBN.IBNcellID));
end

for i = 1:size(SaccadicToIntegrator.SaccadicToSaccadicABD,2)
    index = find(SaccadicToIntegrator.SaccadicToSaccadicABD(i) == SaccadicAxons);
    SaccadicToIBN.ABDcounts(i,:) = histcounts(categorical(Sacc(index).Outputs),categorical(SaccadicToIBN.IBNcellID));
end

for i = 1:size(SaccadicToIntegrator.SaccadicToSaccadicABDi,2)
    index = find(SaccadicToIntegrator.SaccadicToSaccadicABDi(i) == SaccadicAxons);
    SaccadicToIBN.ABDicounts(i,:) = histcounts(categorical(Sacc(index).Outputs),categorical(SaccadicToIBN.IBNcellID));
end

 SaccadicToIBN.AllABDcounts = sum(SaccadicToIBN.ABDcounts,1);
 SaccadicToIBN.AllABDicounts = sum(SaccadicToIBN.ABDicounts,1);
 
 SaccadicToIBN.NormalizedCounts = (SaccadicToIBN.AllABDcounts -  SaccadicToIBN.AllABDicounts) ./ (SaccadicToIBN.AllABDcounts + SaccadicToIBN.AllABDicounts)

 figure;
 subplot(4,4,1)
 histogram(SaccadicToIBN.NormalizedCounts,-1:0.1:1,'FaceColor','k');
 box off;
 axis square;
 xlabel('(Sa-Sai)/(Sa+Sai)');
 
 figure;
 SaccadicToIBN.ABDcellIDs = SaccadicToIBN.IBNcellID(SaccadicToIBN.NormalizedCounts>0);
 SaccadicToIBN.ABDicellIDs = SaccadicToIBN.IBNcellID(SaccadicToIBN.NormalizedCounts<0);
 
 figure;
 subplot(1,2,1)
 transform_swc_AV(SaccadicToIBN.ABDcellIDs,IbnABDcolor,[],true,true);
 subplot(1,2,2)
 transform_swc_AV(SaccadicToIBN.ABDicellIDs,IbnABDicolor,[],true,true);



%% motor patterns of integrator neurons

IntegratorsSaccABD = SaccadicToIntegrator.SaccABDIntCellIDs;
IntegratorsSaccABDi = SaccadicToIntegrator.SaccABDiIntCellIDs;

SaccadicToIntegrator.SaccABDIntMotorDist = isMotor(SaccadicToIntegrator.SaccABDIntCellIDs,df);
SaccadicToIntegrator.SaccABDiIntMotorDist = isMotor(SaccadicToIntegrator.SaccABDiIntCellIDs,df);

figure;
subplot(4,4,1)
imagesc(SaccadicToIntegrator.SaccABDIntMotorDist(:,2:end));
subplot(4,4,2)
imagesc(SaccadicToIntegrator.SaccABDiIntMotorDist(:,2:end));

for i = 1:numel(IntegratorsSaccABD)
    Int(i) = InputsByClass(IntegratorsSaccABD(i),df);
    
    for j = 1:numel(Int(i).Saccadic)
    end
end


%% Saccade --> Saccade interactions

% unbiased saccadic counts
for i = 1:size(SaccadicAxons,2)
    SaccadicToIntegrator.SaccadicCount(i,:) = histcounts(categorical(Sacc(i).Saccadic),categorical(SaccadicToIntegrator.SaccadicCellID));
end

SaccadicToIntegrator.SaccadicCount = SaccadicToIntegrator.SaccadicCount.*~eye(size(SaccadicToIntegrator.SaccadicCount));
SaccadicToIntegrator.SaccadicToSaccadicABD = SaccadicToIntegrator.SaccadicCellID(SaccadicToIntegrator.ABDorder);
SaccadicToIntegrator.SaccadicToSaccadicABDi = SaccadicToIntegrator.SaccadicCellID(SaccadicToIntegrator.ABDiorder);

figure;
imagesc(SaccadicToIntegrator.SaccadicCount([SaccadicToIntegrator.ABDorder;SaccadicToIntegrator.ABDiorder],[SaccadicToIntegrator.ABDorder;SaccadicToIntegrator.ABDiorder]))
colormap(colorcet('L8'));
line([0,size(SaccadicToIntegrator.SaccadicCount,1)],[size(SaccadicToIntegrator.ABDorder,1),size(SaccadicToIntegrator.ABDorder,1)],'color','w');
line([size(SaccadicToIntegrator.ABDorder,1),size(SaccadicToIntegrator.ABDorder,1)],[0,size(SaccadicToIntegrator.SaccadicCount,1)],'color','w');

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

vestibularCellIds = [76631,76700,76656,76655,76660,76679,77393,77395,77375,...
     77815,77803,77807,80591,80579,77255,77697,77364,80223,81126,80760,81117,...
     80672,80622,80572,80569,81122,78426,81328,81575,81582,81870,82161,82165,81553,82169]; 
 
 MVNs = [76664 76784 76783 77394 77326 77772 77766 77756 77753 77767 77775 ...
     77780 80883 80247 80669 80698 80762 80748 80696 80761 80749 81070 81044 ...
     80991 80967 81008 80883 81045 80986 81023 78647 78619 80247 81141];
 
 lateralVSaccadic = [80163 80167 80177 80179 76688 80204 80206 80210 76682];
 
bushySaccadicLateral = [79067 79058 79069 79072 79055 79074 79062 79077 79085 ...
    79080 79086 79064 78576 78540 77146 78646 78542 78541 78543 80728]; 


%sort Saccadic

sortedSaccadics = SaccadicMotorDistribution.normalized(flip(sortSac),1);


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
ABDRestPop.MotorDist = isMotor(ABDRestPop.cellIDs,df);
ABDRestPop.OSI = (sum(ABDRestPop.MotorDist(:,2:3),2)-sum(ABDRestPop.MotorDist(:,4:5),2)) ./ ...
    (sum(ABDRestPop.MotorDist(:,2:3),2)+sum(ABDRestPop.MotorDist(:,4:5),2));
[~,ABDRestPop.OSIsortOrder] = sort(ABDRestPop.OSI);
ABDRestPop.cellIDOSIOrdered = ABDRestPop.cellIDs(flip(ABDRestPop.OSIsortOrder));

 %%
%Integrator, Saccadic, vestibular, Abducens
%  MatOrder = [ SaccadicToIntegrator.SaccadicToSaccadicABD'; SaccadicToIntegrator.SaccadicToSaccadicABDi'; ...
%               SaccadicToIBN.ABDcellIDs';SaccadicToIBN.ABDicellIDs';...
%               vestibularCellIds'; MVNs';...
%              SaccadicToIntegrator.SaccABDIntCellIDs;SaccadicToIntegrator.SaccABDiIntCellIDs;...
%              ABDr_CellIDs';ABDc_CellIDs';ABDIr_CellIDs';ABDIc_CellIDs';...
%              ABDContraPop.cellIDrhombomereOrdered; ABDiContraPop.cellIDrhombomereOrdered ;...
%             ABDRestPop;ABDiRestPop;...
%              ];
%        
% MatOrder = [
%             SaccadicToIntegrator.SaccadicToSaccadicABD'; SaccadicToIntegrator.SaccadicToSaccadicABDi';...
%             %lateralVSaccadic';bushySaccadicLateral';...
%             vestibularCellIds'; MVNs';...
%             SaccadicToIntegrator.SaccABDIntCellIDs; SaccadicToIntegrator.SaccABDiIntCellIDs;...
%             SaccadicToIBN.ABDcellIDs';SaccadicToIBN.ABDicellIDs';...
%             ABDr_CellIDs';ABDc_CellIDs';ABDIr_CellIDs';ABDIc_CellIDs'];
%        

MatOrder = [ sortedSaccadics;
    %SaccadicToIntegrator.SaccadicToSaccadicABD';SaccadicToIntegrator.SaccadicToSaccadicABDi';...
      lateralVSaccadic';bushySaccadicLateral';...
    SaccadicToIBN.ABDcellIDs';SaccadicToIBN.ABDicellIDs';...
    vestibularCellIds'; MVNs';...
    confirmedALX';putativeALX';confirmedDBX';putativeDBX';...
    ABDr_CellIDs';ABDc_CellIDs';ABDIr_CellIDs';ABDIc_CellIDs';...
    ABDContraPop.cellIDOSIrhombomereOrdered;...
    ABDRestPop.cellIDOSIOrdered;...
    ];

        
        
for i = 1:numel(MatOrder)
    MatIndex(i) = find(AllCells == MatOrder(i),1);
end

%subplot(2,2,1)
connMat = ConnMatrixPre(MatIndex,MatIndex);
%cspy(connMat,'Colormap',colorcet('L2'),'Levels',255,'MarkerSize',10);
%subplot(2,2,1);
%imagesc(connMat);
cspy(connMat,'Colormap',colorcet('R3','N',15),'Levels',15,'MarkerSize',7);
axis square;
%colormap(colorcet('L17','N',15,'reverse',1));
lighting phong;
material shiny;
set(gca, 'Xtick',[],'XTickLabel',[],'YTick',[],'YTickLabel',[]);
set(gca,'color',[0.9,0.9,0.9]);
colorbar;
box on;

endABDInt = find(SaccadicToIntegrator.SaccABDIntCellIDs(end) == AllCells(MatIndex))+0.5;
endABDiInt = find(SaccadicToIntegrator.SaccABDiIntCellIDs(end) == AllCells(MatIndex))+0.5;
endSacVent = find(lateralVSaccadic(end) == AllCells(MatIndex))+0.5;
endBushyS = find(bushySaccadicLateral(end) == AllCells(MatIndex))+0.5;
endABDInt = find(putativeALX(end) == AllCells(MatIndex))+0.5;
endABDiInt = find(putativeDBX(end) == AllCells(MatIndex))+0.5;
%endSacABD = find(SaccadicToIntegrator.SaccadicToSaccadicABD(end) == AllCells(MatIndex))+0.5;
endSacABD = find(sortedSaccadics(end) == AllCells(MatIndex))+0.5;
endSACABDi = find(SaccadicToIntegrator.SaccadicToSaccadicABDi(end) == AllCells(MatIndex))+0.5;
endVestABD = find(vestibularCellIds(end) == AllCells(MatIndex))+0.5;
endVestABDi = find(MVNs(end) == AllCells(MatIndex))+0.5;
endIBNordered = find(SaccadicToIBN.ABDcellIDs(end) == AllCells(MatIndex))+0.5;
endIBNrem = find(SaccadicToIBN.ABDicellIDs(end) == AllCells(MatIndex))+0.5;
endContraABD = find(ABDContraPop.cellIDOSIrhombomereOrdered(end) ==  AllCells(MatIndex))+0.5;
%endContraABDi = find( ABDiContraPop.cellIDrhombomereOrdered(end) ==  AllCells(MatIndex))+0.5;
endRestABD = find(ABDRestPop.cellIDOSIOrdered(end) ==  AllCells(MatIndex))+0.5;
%endRestABDi = find(ABDiRestPop(end) ==  AllCells(MatIndex))+0.5;
endABD = find(ABDc_CellIDs(end) == AllCells(MatIndex))+0.5;
endABDi = find(ABDIc_CellIDs(end) == AllCells(MatIndex))+0.5;

hold on;
 line([0,size(connMat,1)],[endSacVent,endSacVent],'color',[0.5,0.5,0.5]);
 line([0,size(connMat,1)],[endBushyS(1),endBushyS(1)],'color',[0.5,0.5,0.5]);
line([0,size(connMat,1)],[endABDInt(1),endABDInt(1)],'color',[0.5,0.5,0.5]);
line([0,size(connMat,1)],[endABDiInt,endABDiInt],'color',[0.5,0.5,0.5]);
line([0,size(connMat,1)],[endSacABD,endSacABD],'color',[0.5,0.5,0.5]);
%line([0,size(connMat,1)],[endSACABDi,endSACABDi],'color',[0.5,0.5,0.5]);
 line([0,size(connMat,1)],[endVestABD(1),endVestABD(1)],'color',[0.5,0.5,0.5]);
 line([0,size(connMat,1)],[endVestABDi,endVestABDi],'color',[0.5,0.5,0.5]);
line([0,size(connMat,1)],[endIBNordered,endIBNordered],'color',[0.5,0.5,0.5]);
line([0,size(connMat,1)],[endIBNrem,endIBNrem],'color',[0.5,0.5,0.5]);
  line([0,size(connMat,1)],[endContraABD,endContraABD],'color','k');
%  line([0,size(connMat,1)],[endContraABDi,endContraABDi],'color','k');
  line([0,size(connMat,1)],[endRestABD,endRestABD],'color','k');
%  line([0,size(connMat,1)],[endRestABDi(2),endRestABDi(2)],'color','k');
line([0,size(connMat,1)],[endABD,endABD],'color',[0.5,0.5,0.5]);
line([0,size(connMat,1)],[endABDi,endABDi],'color',[0.5,0.5,0.5]);


% 
 line([endSacVent,endSacVent],[0,size(connMat,1)],'color',[0.5,0.5,0.5]);
 line([endBushyS(1),endBushyS(1)],[0,size(connMat,1)],'color',[0.5,0.5,0.5]);
  line([endABDInt(1),endABDInt(1)],[0,size(connMat,1)],'color',[0.5,0.5,0.5]);
  line([endABDiInt,endABDiInt],[0,size(connMat,1)],'color',[0.5,0.5,0.5]);
line([endSacABD,endSacABD],[0,size(connMat,1)],'color',[0.5,0.5,0.5]);
%line([endSACABDi,endSACABDi],[0,size(connMat,1)],'color',[0.5,0.5,0.5]);
line([endVestABD(1),endVestABD(1)],[0,size(connMat,1)],'color',[0.5,0.5,0.5]);
line([endVestABDi,endVestABDi],[0,size(connMat,1)],'color',[0.5,0.5,0.5]);
line([endIBNordered,endIBNordered],[0,size(connMat,1)],'color',[0.5,0.5,0.5]);
line([endIBNrem,endIBNrem],[0,size(connMat,1)],'color',[0.5,0.5,0.5]);
 line([endContraABD,endContraABD],[0,size(connMat,1)],'color','k');
%  line([endContraABDi,endContraABDi],[0,size(connMat,1)],'color','k');
 line([endRestABD,endRestABD],[0,size(connMat,1)],'color','k');
%  line([endRestABDi(2),endRestABDi(2)],[0,size(connMat,1)],'color','k');
line([endABD,endABD],[0,size(connMat,1)],'color',[0.5,0.5,0.5]);
line([endABDi,endABDi],[0,size(connMat,1)],'color',[0.5,0.5,0.5]);

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
    bushySaccadicLateralSm(i) = sum(ismember(a,SaccadicProjectingToABDexclusively));
    bushySaccadicLateralSi(i) = sum(ismember(a,SaccadicProjectingToABDiexclusively));
    clear a;
end
