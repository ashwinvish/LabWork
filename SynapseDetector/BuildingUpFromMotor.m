clear

df = readtable('/Users/ashwin/Documents/SynapseDetector/04152019.csv');

% load Saccadic neurons

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

allSaccade = [triangle,sparseSaccadic,bushySaccadicMedial,bushySaccadicLateral,...
    putativeBushySaccadic,lateralDSaccadic,lateralVSaccadic,unknownSaccadic];

% load Abducens neurons


ABDr_CellIDs = [77648 77710 77300 77705 77305 77301 77709 77672 77302 82194 82192 82193 82146 82145 82143 82140];
ABDc_CellIDs = [77154 77646 77682 77628 77295 77652 77292 77688 77654 77658 77657 77662 77296 81172 82195 82196 82212 82213 82197];
ABDIr_CellIDs = [77631, 77150, 77618, 77886, 78547, 77158, 78556, 78552, 77665, 77668, 77634, 78553];
ABDIc_CellIDs = [77148, 77625, 77641, 77692, 77144, 77643, 77640, 79051, 79066, 78574];

allMotor = [ABDr_CellIDs,ABDc_CellIDs];
allInterNuclear = [ABDIr_CellIDs,ABDIc_CellIDs];
allABD = [allMotor,allInterNuclear];

load('ABDVols.mat') ; % load volume measures for Abducens neurons

%% 
for i = 1:numel(allMotor)
    ABD(i) = InputsByClass(allMotor(i),df);
end

for i = 1:numel(allMotor)
    temp = find(ABDvols(:,1)==allMotor(i));
    if ~isempty(temp)
        ABD(i).volume = ABDvols(temp,2);
    else
        ABD(i).volume = [];
    end
end

for i = 1:numel(allInterNuclear)
        ABDi(i) = InputsByClass(allInterNuclear(i),df);
end

for i = 1:numel(allInterNuclear)
    temp = find(ABDvols(:,1)==allInterNuclear(i));
    if ~isempty(temp)
        ABDi(i).volume = ABDvols(temp,2);
    else
        ABDi(i).volume = [];
    end
end

%% 

a = find(ismember(allMotor,ABDr_CellIDs));
b = find(ismember(allMotor,ABDc_CellIDs));

for i = 1:numel(a)
    scatter3(ABD(a(i)).Origin(1),ABD(a(i)).Origin(3),ABD(a(i)).volume ./1e9);
    hold on;
end

for i = 1:numel(b)
    scatter3(ABD(b(i)).Origin(1),ABD(b(i)).Origin(3),ABD(b(i)).volume ./1e9);
    hold on;
end

%%
[ABDsort,ABDvolOrder] = sort([ABD.volume]);
ABDvolOrder(isnan(ABDsort)) = [];
[~,ABDivolOrder] = sort([ABDi.volume]);

for i = 1:numel(ABD)
    ABDSynapseNumbers(i) = size(ABD(i).Saccadic,1);
end

plot(1:length(ABDSynapseNumbers), ABDSynapseNumbers,'ko');
showfit(ezfit(1:length(ABDSynapseNumbers), ABDSynapseNumbers,'affine'),'fitcolor','k');
hold on;
plot(1:length(ABDSynapseNumbers(ABDvolOrder)), ABDSynapseNumbers(ABDvolOrder),'ro');
showfit(ezfit(1:length(ABDSynapseNumbers(ABDvolOrder)), ABDSynapseNumbers(ABDvolOrder),'affine'),'fitcolor','r');

for i = 1:numel(ABDSynapseNumbers(ABDvolOrder))
    if ~isempty(ABD(ABDvolOrder(i)).Tree)
    subplot(6,6,i)
    histfit(ABD(ABDvolOrder(i)).PathLength(ABD(ABDvolOrder(i)).isContra)./max(Pvec_tree(ABD(ABDvolOrder(i)).Tree{1})),10,'kernel')
    end
end

ABDsorterOrigin = vertcat(ABD(ABDvolOrder).Origin);
ABDsortedVolume = [ABD(ABDvolOrder).volume];
ABDsortedVolume(isnan(ABDsortedVolume)) = [];
ABDsortedCellID = [ABD(ABDvolOrder).cellID]'

scatter3(ABDsorterOrigin(:,1),ABDsorterOrigin(:,3),ABDsortedVolume'./1e9,20,'filled');
text(ABDsorterOrigin(:,1),ABDsorterOrigin(:,3),ABDsortedVolume'./1e9, num2str(ABDsortedCellID));
%%

for i = 1:numel(allSaccade)
    motorDiff(i,:) = isMotor(allSaccade(i),df);
end

ABDheavy = allSaccade(find((motorDiff(:,2)+motorDiff(:,3)) > (motorDiff(:,4)+motorDiff(:,5))));
ABDiheavy = allSaccade(find((motorDiff(:,2)+motorDiff(:,3)) < (motorDiff(:,4)+motorDiff(:,5))));

%%

colors = colorcet('D2','N',5);
subplot(2,3,[1,4])
transform_swc_AV(ABDheavy,colors(1,:),[],true,false);
transform_swc_AV(ABDiheavy,colors(5,:),[],false,false);

ABDHeavyOrigins = getOrigin(ABDheavy);
ABDiHeavyOrigins = getOrigin(ABDiheavy);


ABDHeavyFit = fitdist(ABDHeavyOrigins(:,2),'Kernel','BandWidth',4);
ABDx = 0:2:1280;
ABDy = pdf(ABDHeavyFit,ABDx);

ABDiHeavyFit = fitdist(ABDiHeavyOrigins(:,2),'Kernel','BandWidth',4);
ABDix = 0:2:1280;
ABDiy = pdf(ABDiHeavyFit,ABDix);

%
% subplot(2,3,[2,5])
%
% %histogram(ABDHeavyOrigins(:,2),Bins,'EdgeColor',colors(1,:),'Orientation','horizontal','DisplayStyle','stairs','LineWidth',2);
% plot(ABDx,ABDy,'color',colors(1,:),'LineWidth',2);
% hold on;
% %histogram(ABDiHeavyOrigins(:,2),Bins,'EdgeColor',colors(5,:),'Orientation','horizontal','DisplayStyle','stairs','LineWidth',2);
%
% %set(gca,'YDir','reverse','YLim',[0,1280]);
% plot(ABDix,ABDiy,'color',colors(5,:),'LineWidth',2);
%
% daspect([1,30,1]);


