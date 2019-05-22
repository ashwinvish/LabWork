%clear;

colorPallete = cbrewer('div','BrBG',5);

smallColor = colorPallete(1,:);
largeColor = colorPallete(5,:);

small = 'MIF';
large = 'SIF';

ABDr_CellIDs = [77648 77710 77300 77705 77305 77301 77709 77672 77302 82194 82192 82193 82146 82145 82143 82140];
ABDc_CellIDs = [77154 77646 77682 77628 77295 77652 77292 77688 77654 77658 77657 77662 77296 81172 82195 82196 82212 82213 82197];
ABDIr_CellIDs = [77631, 77150, 77618, 77886, 78547, 77158, 78556, 78552, 77665, 77668, 77634, 78553];
ABDIc_CellIDs = [77148, 77625, 77641, 77692, 77144, 77643, 77640, 79051, 79066, 78574];

AllABD = [ABDr_CellIDs,ABDc_CellIDs,ABDIr_CellIDs,ABDIc_CellIDs];
Allmotor = [ABDr_CellIDs,ABDc_CellIDs];

% AllABDradiusIndex = [2,
% AllABDradius = [2400,];
load('ABDVols.mat');
for i = 1:numel(AllABD)
    if ismember(AllABD(i),ABDvols(:,1))
        index = find(AllABD(i)==ABDvols(:,1));
        ABD.vol(i) = ABDvols(index,2);
        ABD.cellID(i) = AllABD(i);
    end
end

vol = vol./1e9;% convert to microns
vol(vol==0) = NaN;

figure;
subplot(4,4,1)
histogram(vol(1:numel(Allmotor)),20,'FaceColor','g');
title('ABD volumes');
box off;
line([110,110],[0,5],'color','k','LineStyle','--','LineWidth',2);
axis square;

subplot(4,4,2)
histogram(vol(1+numel(Allmotor):end),20,'FaceColor','m');
title('ABDi volumes');
box off;
axis square;

load('ABDr.mat');
load('ABDc.mat');
load('ABDIr.mat');
load('ABDIc.mat');


for i = 1:numel(AllABD)
    if ismember (AllABD(i),ABDr_CellIDs)
        index = find(AllABD(i)==ABDr_CellIDs);
        ABDOrigin(i,:) = ABDr(index).Origin;
    elseif ismember (AllABD(i),ABDc_CellIDs)
        index = find(AllABD(i)==ABDc_CellIDs);
        ABDOrigin(i,:) = ABDc(index).Origin;
    elseif ismember (AllABD(i),ABDIr_CellIDs)
        index = find(AllABD(i)==ABDIr_CellIDs);
        ABDOrigin(i,:) = ABDIr(index).Origin;
    else
        index = find(AllABD(i)==ABDIc_CellIDs);
        ABDOrigin(i,:) = ABDIc(index).Origin;
    end
end


subplot(4,4,3)

scatter3(ABDOrigin(:,1),ABDOrigin(:,3),vol',20,'o','MarkerFaceColor',smallColor,'MarkerEdgeColor','none');

xlabel('M<-->L (\mu)');
ylabel('D<-->V (\mu)');
zlabel('vol (\mu^3)');
set(gca,'YDir','reverse');

axis square


% subplot(4,4,5)
% histogram(vol2(1:22),10,'FaceColor','g');
% title('ABD volumes');
% box off;
% line([110,110],[0,10],'color','k','LineStyle','--','LineWidth',2);
% axis square;
%
% subplot(4,4,6)
% histogram(vol2(23:end),10,'FaceColor','m');
% title('ABDi volumes');
% box off;
% axis square;


prompt = 'Plot in correct order? Y/N [Y]:';
str = input(prompt,'s');
if isempty(str)
    str = 'Y';
end


%% plot in correct order
if (str == 'Y' )
    
    
    
    % correct order
    %smallABDCellIDs = Allmotor(find(vol(1:22)<110));
    smallABDCellIDs = [82140,82145,82143,77648,82146,77302,82213,77654];
    save('smallABDneurons.mat','smallABDCellIDs');
    %smallABDvol = vol(find(vol(1:22)<132));
    %largeABDCellIDs = Allmotor(find(vol(1:22)>110));
    largeABDCellIDs= [77301, 82194,77705,77710,77305,77672,77300,77709, 77646,82212,77295,81172,77657,82197,77682,77154, 77652,77658,77628,77292,77688,82195,77296];
    save('largeABDneurons.mat','largeABDCellIDs');
    %largeABDvol = vol(find(vol(1:22)>132));
    
    
    for i = 1:numel(ABDr_CellIDs)
        if ismember(ABDr_CellIDs(i),smallABDCellIDs)
            smallABD(i).cellID = ABDr_CellIDs(i);
            smallABD(i).stats = ABDr(i);
            
        elseif ismember(ABDr_CellIDs(i),largeABDCellIDs);
            largeABD(i).cellID = ABDr_CellIDs(i);
            largeABD(i).stats = ABDr(i);
        end
    end
    
    j = numel(ABDr_CellIDs)+1;
    
    for i = 1:numel(ABDc_CellIDs)
        if ismember(ABDc_CellIDs(i),smallABDCellIDs)
            smallABD(i+j).cellID = ABDc_CellIDs(i);
            smallABD(i+j).stats = ABDc(i);
        elseif ismember(ABDc_CellIDs(i),largeABDCellIDs);
            largeABD(i+j).cellID = ABDc_CellIDs(i);
            largeABD(i+j).stats = ABDc(i);
        end
    end
    
    
    index = 1;
    smallSaccadicAxons = [];
    smallContraAxons = [];
    smallVestibularAxons = [];
    smallIntegratorAxons = [];
    smallABDorigin = [];
    
    for i = 1:numel(smallABD)
        if ~isempty(smallABD(i).cellID)
            numberOfSmallSaccadicSynapses(index) = size(smallABD(i).stats.Saccadic,1);
            numberOfSmallVestibularSynapses(index) = size(smallABD(i).stats.Vestibular,1);
            smallSaccadicAxons = [smallSaccadicAxons;smallABD(i).stats.Saccadic];
            smallContraAxons = [smallContraAxons;smallABD(i).stats.Contra];
            smallVestibularAxons = [smallVestibularAxons;smallABD(i).stats.Vestibular];
            smallIntegratorAxons = [smallIntegratorAxons;smallABD(i).stats.Integrator];
            smallABDorigin = [smallABDorigin;smallABD(i).stats.Origin];
            smallSynapsesPerABD(1:size(unique(smallABD(i).stats.Saccadic),1),index) = histc(smallABD(i).stats.Saccadic,unique(smallABD(i).stats.Saccadic));
            
            smallABD_SaccadicDist(index,:) = histcounts(smallABD(i).stats.SaccadicDist,0:0.1:1,'Normalization','probability');
            smallABD_ContraDist(index,:) = histcounts(smallABD(i).stats.ContraDist,0:0.1:1,'Normalization','probability');
            smallABD_VestibularDist(index,:) = histcounts(smallABD(i).stats.VestibularDist,0:0.1:1,'Normalization','probability');
            smallABD_IntegratorDist(index,:) = histcounts(smallABD(i).stats.IntegratorDist,0:0.1:1,'Normalization','probability');
            smallABD_EverythingElseDist(index,:) = histcounts(smallABD(i).stats.EverythingElseDist,0:0.1:1,'Normalization','probability');
            index = index+1;
        end
    end
    
    index = 1;
    largeSaccadicAxons = [];
    largeVestibularAxons = [];
    largeContraAxons = [];
    largeIntegratorAxons = [];
    largeABDorigin = [];
    
    for i = 1:numel(largeABD)
        if ~isempty(largeABD(i).cellID)
            numberOfLargeSaccadicSynapses(index) = size(largeABD(i).stats.Saccadic,1);
            numberOfLargeVestibularSynapses(index) = size(largeABD(i).stats.Vestibular,1);
            largeSaccadicAxons = [largeSaccadicAxons;largeABD(i).stats.Saccadic];
            largeVestibularAxons =[largeVestibularAxons;largeABD(i).stats.Vestibular];
            largeContraAxons = [largeContraAxons;largeABD(i).stats.Contra];
            largeIntegratorAxons =[largeIntegratorAxons;largeABD(i).stats.Integrator];
            largeABDorigin = [largeABDorigin ; largeABD(i).stats.Origin];
            largeSynapsesPerABD(1:size(unique(largeABD(i).stats.Saccadic),1),index) = histc(largeABD(i).stats.Saccadic,unique(largeABD(i).stats.Saccadic));
            
            largeABD_SaccadicDist(index,:) = histcounts(largeABD(i).stats.SaccadicDist,0:0.1:1,'Normalization','probability');
            largeABD_ContraDist(index,:) = histcounts(largeABD(i).stats.ContraDist,0:0.1:1,'Normalization','probability');
            largeABD_VestibularDist(index,:) = histcounts(largeABD(i).stats.VestibularDist,0:0.1:1,'Normalization','probability');
            largeABD_IntegratorDist(index,:) = histcounts(largeABD(i).stats.IntegratorDist,0:0.1:1,'Normalization','probability');
            largeABD_EverythingElseDist(index,:) = histcounts(largeABD(i).stats.EverythingElseDist,0:0.1:1,'Normalization','probability');
            index = index+1;
        end
    end
    
    
    %% plot correct order with standard Error
    
    % saccadic dist
    %figure;
    
    
    smallABDnumber = length(smallABD_SaccadicDist);
    largeABDnumber = length(largeABD_SaccadicDist);
    
    % All Inputs
    subplot(4,4,4)
    h1 = shadedErrorBar([0.1:0.1:1],mean([smallABD_SaccadicDist;smallABD_VestibularDist;smallABD_ContraDist;smallABD_IntegratorDist;smallABD_EverythingElseDist])...
        ,std([smallABD_SaccadicDist;smallABD_VestibularDist;smallABD_ContraDist;smallABD_IntegratorDist;smallABD_EverythingElseDist])./sqrt(smallABDnumber),...
        'lineProps',{'Color',smallColor,'LineWidth',2});
    hold on;
    h2 = shadedErrorBar([0.1:0.1:1],mean([largeABD_SaccadicDist;largeABD_VestibularDist;largeABD_ContraDist;largeABD_IntegratorDist;largeABD_EverythingElseDist])...
        ,std([largeABD_SaccadicDist;largeABD_VestibularDist;largeABD_ContraDist;largeABD_IntegratorDist;largeABD_EverythingElseDist])./sqrt(largeABDnumber)...
        ,'lineProps',{'Color',largeColor,'LineWidth',2});
    box off;
    axis square;
    title('All inputs');
    
    legend([h1.patch.Parent], {'Small (mif)','Large (SIF)'},'Location','bestoutside');
    
    % All inputs randomized
    
    % Saccadic
    subplot(4,4,5)
    shadedErrorBar([0.1:0.1:1],mean(smallABD_SaccadicDist),std(smallABD_SaccadicDist)./sqrt(smallABDnumber),'lineProps',{'Color',smallColor,'LineWidth',2});
    hold on;
    shadedErrorBar([0.1:0.1:1],mean(largeABD_SaccadicDist),std(largeABD_SaccadicDist)./sqrt(largeABDnumber),'lineProps',{'Color',largeColor,'LineWidth',2});
    
    box off;
    axis square;
    title('Saccadic');
    
    
    % Vestibular dist
    subplot(4,4,6)
    shadedErrorBar([0.1:0.1:1],mean(smallABD_VestibularDist),std(smallABD_VestibularDist)./sqrt(smallABDnumber),'lineProps',{'Color',smallColor,'LineWidth',2});
    hold on;
    shadedErrorBar([0.1:0.1:1],mean(largeABD_VestibularDist),std(largeABD_VestibularDist)./sqrt(largeABDnumber),'lineProps',{'Color',largeColor,'LineWidth',2});
    
    box off;
    axis square;
    title('Vestibular');
    
    % contra Dist
    subplot(4,4,7)
    shadedErrorBar([0.1:0.1:1],mean(smallABD_ContraDist),std(smallABD_ContraDist)./sqrt(smallABDnumber),'lineProps',{'Color',smallColor,'LineWidth',2});
    hold on;
    shadedErrorBar([0.1:0.1:1],mean(largeABD_ContraDist),std(largeABD_ContraDist)./sqrt(largeABDnumber),'lineProps',{'Color',largeColor,'LineWidth',2});
    
    box off;
    axis square;
    title('Contra');
    
    % integrator Dist
    subplot(4,4,8)
    shadedErrorBar([0.1:0.1:1],mean(smallABD_IntegratorDist),std(smallABD_IntegratorDist)./sqrt(smallABDnumber),'lineProps',{'Color',smallColor,'LineWidth',2});
    hold on;
    shadedErrorBar([0.1:0.1:1],mean(largeABD_IntegratorDist),std(largeABD_IntegratorDist)./sqrt(largeABDnumber),'lineProps',{'Color',largeColor,'LineWidth',2});
    
    box off;
    axis square;
    title('Integrator');
    
else
    
    % plot with randomized order
    
    load('ABDr.mat');
    load('ABDc.mat');
    load('ABDIr.mat');
    load('ABDIc.mat');
    
    % randomized order
    
    smallABDCellIDs = Allmotor(randi(22,[1,11]));
    largeABDCellIDs = Allmotor(randi(22,[1,11]));
    
    for i = 1:numel(ABDr_CellIDs)
        if ismember(ABDr_CellIDs(i),smallABDCellIDs)
            smallABD(i).cellID = ABDr_CellIDs(i);
            smallABD(i).stats = ABDr(i);
            
        elseif ismember(ABDr_CellIDs(i),largeABDCellIDs);
            largeABD(i).cellID = ABDr_CellIDs(i);
            largeABD(i).stats = ABDr(i);
        end
    end
    
    j = numel(ABDr_CellIDs)+1;
    
    for i = 1:numel(ABDc_CellIDs)
        if ismember(ABDc_CellIDs(i),smallABDCellIDs)
            smallABD(i+j).cellID = ABDc_CellIDs(i);
            smallABD(i+j).stats = ABDc(i);
        elseif ismember(ABDc_CellIDs(i),largeABDCellIDs);
            largeABD(i+j).cellID = ABDc_CellIDs(i);
            largeABD(i+j).stats = ABDc(i);
        end
    end
    
    index = 1;
    smallSaccadicAxons = [];
    smallContraAxons = [];
    smallVestibularAxons = [];
    smallIntegratorAxons = [];
    
    for i = 1:numel(smallABD)
        if ~isempty(smallABD(i).cellID)
            numberOfSmallSaccadicSynapses(index) = size(smallABD(i).stats.Saccadic,1);
            numberOfSmallVestibularSynapses(index) = size(smallABD(i).stats.Vestibular,1);
            smallSaccadicAxons = [smallSaccadicAxons;smallABD(i).stats.Saccadic];
            smallContraAxons = [smallContraAxons;smallABD(i).stats.Contra];
            smallVestibularAxons = [smallVestibularAxons;smallABD(i).stats.Vestibular];
            smallIntegratorAxons = [smallIntegratorAxons;smallABD(i).stats.Integrator];
            
            smallSynapsesPerABD(1:size(unique(smallABD(i).stats.Saccadic),1),index) = histc(smallABD(i).stats.Saccadic,unique(smallABD(i).stats.Saccadic));
            
            smallABD_SaccadicDist(index,:) = histcounts(smallABD(i).stats.SaccadicDist,0:0.1:1,'Normalization','probability');
            smallABD_ContraDist(index,:) = histcounts(smallABD(i).stats.ContraDist,0:0.1:1,'Normalization','probability');
            smallABD_VestibularDist(index,:) = histcounts(smallABD(i).stats.VestibularDist,0:0.1:1,'Normalization','probability');
            smallABD_IntegratorDist(index,:) = histcounts(smallABD(i).stats.IntegratorDist,0:0.1:1,'Normalization','probability');
            smallABD_EverythingElseDist(index,:) = histcounts(smallABD(i).stats.EverythingElseDist,0:0.1:1,'Normalization','probability');
            index = index+1;
        end
    end
    
    
    index = 1;
    largeSaccadicAxons = [];
    largeVestibularAxons = [];
    largeContraAxons = [];
    largeIntegratorAxons = [];
    
    for i = 1:numel(largeABD)
        if ~isempty(largeABD(i).cellID)
            numberOfLargeSaccadicSynapses(index) = size(largeABD(i).stats.Saccadic,1);
            numberOfLargeVestibularSynapses(index) = size(largeABD(i).stats.Vestibular,1);
            largeSaccadicAxons = [largeSaccadicAxons;largeABD(i).stats.Saccadic];
            largeVestibularAxons =[largeVestibularAxons;largeABD(i).stats.Vestibular];
            largeContraAxons = [largeContraAxons;largeABD(i).stats.Contra];
            largeIntegratorAxons =[largeIntegratorAxons;largeABD(i).stats.Integrator];
            
            largeSynapsesPerABD(1:size(unique(largeABD(i).stats.Saccadic),1),index) = histc(largeABD(i).stats.Saccadic,unique(largeABD(i).stats.Saccadic));
            
            largeABD_SaccadicDist(index,:) = histcounts(largeABD(i).stats.SaccadicDist,0:0.1:1,'Normalization','probability');
            largeABD_ContraDist(index,:) = histcounts(largeABD(i).stats.ContraDist,0:0.1:1,'Normalization','probability');
            largeABD_VestibularDist(index,:) = histcounts(largeABD(i).stats.VestibularDist,0:0.1:1,'Normalization','probability');
            largeABD_IntegratorDist(index,:) = histcounts(largeABD(i).stats.IntegratorDist,0:0.1:1,'Normalization','probability');
            largeABD_EverythingElseDist(index,:) = histcounts(largeABD(i).stats.EverythingElseDist,0:0.1:1,'Normalization','probability');
            index = index+1;
        end
    end
    
    %  plot randomized order with standard Error
    
    % saccadic dist
    %figure;
    colorPallete = cbrewer('div','BrBG',5);
    
    smallColor = colorPallete(1,:);
    largeColor = colorPallete(5,:);
    
    small = 'MIF';
    large = 'SIF';
    
    smallABDnumber = length(smallABD_SaccadicDist);
    largeABDnumber = length(largeABD_SaccadicDist);
    
    % All Inputs
    subplot(4,4,11)
    shadedErrorBar([0.1:0.1:1],mean([smallABD_SaccadicDist;smallABD_VestibularDist;smallABD_ContraDist;smallABD_IntegratorDist;smallABD_EverythingElseDist])...
        ,std([smallABD_SaccadicDist;smallABD_VestibularDist;smallABD_ContraDist;smallABD_IntegratorDist;smallABD_EverythingElseDist])./sqrt(smallABDnumber),...
        'lineProps',{'Color',smallColor,'LineWidth',2});
    hold on;
    shadedErrorBar([0.1:0.1:1],mean([largeABD_SaccadicDist;largeABD_VestibularDist;largeABD_ContraDist;largeABD_IntegratorDist;largeABD_EverythingElseDist])...
        ,std([largeABD_SaccadicDist;largeABD_VestibularDist;largeABD_ContraDist;largeABD_IntegratorDist;largeABD_EverythingElseDist])./sqrt(largeABDnumber)...
        ,'lineProps',{'Color',largeColor,'LineWidth',2});
    box off;
    axis square;
    title('All inputs');
    
    % All inputs randomized
    
    % Saccadic
    subplot(4,4,13)
    shadedErrorBar([0.1:0.1:1],mean(smallABD_SaccadicDist),std(smallABD_SaccadicDist)./sqrt(smallABDnumber),'lineProps',{'Color',smallColor,'LineWidth',2});
    hold on;
    shadedErrorBar([0.1:0.1:1],mean(largeABD_SaccadicDist),std(largeABD_SaccadicDist)./sqrt(largeABDnumber),'lineProps',{'Color',largeColor,'LineWidth',2});
    
    box off;
    axis square;
    title('Saccadic');
    
    
    % Vestibular dist
    subplot(4,4,14)
    shadedErrorBar([0.1:0.1:1],mean(smallABD_VestibularDist),std(smallABD_VestibularDist)./sqrt(smallABDnumber),'lineProps',{'Color',smallColor,'LineWidth',2});
    hold on;
    shadedErrorBar([0.1:0.1:1],mean(largeABD_VestibularDist),std(largeABD_VestibularDist)./sqrt(largeABDnumber),'lineProps',{'Color',largeColor,'LineWidth',2});
    
    box off;
    axis square;
    title('Vestibular');
    
    % contra Dist
    subplot(4,4,15)
    shadedErrorBar([0.1:0.1:1],mean(smallABD_ContraDist),std(smallABD_ContraDist)./sqrt(smallABDnumber),'lineProps',{'Color',smallColor,'LineWidth',2});
    hold on;
    shadedErrorBar([0.1:0.1:1],mean(largeABD_ContraDist),std(largeABD_ContraDist)./sqrt(largeABDnumber),'lineProps',{'Color',largeColor,'LineWidth',2});
    
    box off;
    axis square;
    title('Contra');
    
    % integrator Dist
    subplot(4,4,16)
    shadedErrorBar([0.1:0.1:1],mean(smallABD_IntegratorDist),std(smallABD_IntegratorDist)./sqrt(smallABDnumber),'lineProps',{'Color',smallColor,'LineWidth',2});
    hold on;
    shadedErrorBar([0.1:0.1:1],mean(largeABD_IntegratorDist),std(largeABD_IntegratorDist)./sqrt(largeABDnumber),'lineProps',{'Color',largeColor,'LineWidth',2});
    
    box off;
    axis square;
    title('Integrator');
    
end
% %% plot neurons
%
% % correct order
% smallABDCellIDs = Allmotor(find(vol(1:22)<132));
% largeABDCellIDs = Allmotor(find(vol(1:22)>132));
%
%
% figure;
% transform_swc_AV(smallABDCellIDs,smallColor,[],true,false);
% transform_swc_AV(largeABDCellIDs,largeColor,[],false,false);
%
% figure;
%
% smallOnlySaccadicAxons = setdiff(smallSaccadicAxons,largeSaccadicAxons);
% largeOnlySaccdicAxons = setdiff(largeSaccadicAxons,smallSaccadicAxons);
% commonSaccadicAxons = intersect(smallSaccadicAxons,largeSaccadicAxons);
%
%
% transform_swc_AV(smallOnlySaccadicAxons,smallColor,[],true,false);
% transform_swc_AV(largeOnlySaccdicAxons,largeColor,[],false,false);
% %transform_swc_AV(unique(commonSaccadicAxons),[0.8,0.8,0.8],[],false,false);
%
% figure;
%
% commonContraAxons = intersect(smallContraAxons,largeContraAxons);
% smallOnlyContraAxons = setdiff(smallContraAxons,largeContraAxons);
% largeOnlyContraAxons = setdiff(largeContraAxons,smallContraAxons);
%
% transform_swc_AV(smallOnlyContraAxons,smallColor,[],true,false);
% transform_swc_AV(largeOnlyContraAxons,largeColor,[],false,false);
% %transform_swc_AV(unique(commonContraAxons),[0.8,0.8,0.8],[],false,false);
%
% figure;
%
% commonInteratorAxons = intersect(smallIntegratorAxons,largeIntegratorAxons);
% smallOnlyIntegratorAxons = setdiff(largeIntegratorAxons,largeContraAxons);
% largeOnlyIntegratorAxons = setdiff(largeIntegratorAxons,smallIntegratorAxons);
%
% transform_swc_AV(smallOnlyIntegratorAxons,smallColor,[],true,false);
% transform_swc_AV(largeOnlyIntegratorAxons,largeColor,[],false,false);
% %transform_swc_AV(unique(commonInteratorAxons),[0.8,0.8,0.8],[],false,false);
%
% figure;
%
% commonVestibularAxons = intersect(smallVestibularAxons,largeVestibularAxons);
% smallOnlyVestibularAxons = setdiff(smallVestibularAxons,largeVestibularAxons);
% largeOnlyVestibularAxons = setdiff(largeVestibularAxons,smallVestibularAxons);
%
% transform_swc_AV(unique(smallVestibularAxons),smallColor,[],true,false);
% transform_swc_AV(unique(largeVestibularAxons),largeColor,[],false,false);
% %transform_swc_AV(unique(commonVestibularAxons),[0.8,0.8,0.8],[],false,false);
