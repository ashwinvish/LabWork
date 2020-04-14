% plots for clustering data

% gamma 0.4 gives 2 blocks, block1 is RS and block2 is Int.

load block1_gamma038_WithNull.mat % 76749,79560 , gamma = 0.38
load block2_gamma038_WithNull.mat % 79059,77806, gamma = 0.38

block1_gamma038_WithNull.color = [1,0.5,0];
block2_gamma038_WithNull.color = [0,0.5,1];


% gamma 0.8 gives 6 blocks, block5 is ABDi and block6 is ABDm projecting

% load block1_gamma09.mat % With null model
% load block2_gamma09.mat
% load block3_gamma09.mat
% load block4_gamma09.mat
% load block5_gamma09.mat
% load block6_gamma09.mat

load block1_subMod.mat % onto ABDi
load block2_subMod.mat % onto ABD

load AllCells.mat
load ConnMatrixPre.mat


LoadDataFrame

confirmedALX = [76181 76201 76187 76184 76192 76197];
confirmedDBX = [76199 76200 76182 76183 76185 76186 76189 76191 76188 ];
confirmedBARHL = [76198 76190 76193 76194 76195 76196 ];
confirmedIntegrators = [confirmedALX,confirmedDBX,confirmedBARHL];

temp = distinguishable_colors(10,'w');
block1_subMod_gamm075.color = temp(9,:);
block2_subMod_gamm075.color = temp(8,:);

% %% Add motor neurons
ABDr_CellIDs = [82140,82145,82143,77648,82146,77302,77301,82194,77705,77710,77305,77672,77300,77709,77661];
ABDc_CellIDs = [82213,77654,77646,82212,77295,81172,77657,82197,77682,77154, 77652,77658,77628,77292,77688,82195,77296];
ABDIr_CellIDs = [78553,77150,77886,78547,77631,77158,77665,78552,77668,77618,77634,78556];
ABDIc_CellIDs = [78574,79051,77148,77625,77144,77692,77641,79066,77640,77643];
% 
% % 
% % matOrder_commdet = [block1.cellIDs', block2.cellIDs',...
% %     ABDr_CellIDs,ABDc_CellIDs,ABDIr_CellIDs,ABDIc_CellIDs];
% 
% matOrder_commdet = [AllCellNames(idx)',...
%     ABDr_CellIDs,ABDc_CellIDs,ABDIr_CellIDs,ABDIc_CellIDs];
%    
% [~,matIndex_commdet] = ismember(matOrder_commdet,AllCells);
% gap = 20;
% tempMat = zeros(length(matOrder_commdet)+gap);
% 
% tempMat(1:length(AllCellNames),1:length(AllCellNames)) = ConnMatrixPre(matIndex_commdet(1:length(AllCellNames)),matIndex_commdet(1:length(AllCellNames)));
% tempMat(length(AllCellNames)+gap+1:end,1:length(AllCellNames)) = ConnMatrixPre(matIndex_commdet(length(AllCellNames)+1:end),matIndex_commdet(1:length(AllCellNames)));
% %tempMat(1:length(AllCellNames),length(AllCellNames)+gap+1:end) = ConnMatrixPre(matIndex_commdet(1:length(AllCellNames)),matIndex_commdet(length(AllCellNames)+1:end));
% %tempMat(length(AllCellNames)+gap+1:end,length(AllCellNames)+gap+1:end) = ConnMatrixPre(matIndex_commdet(length(AllCellNames)+1:end),matIndex_commdet(length(AllCellNames)+1:end));
% 
% figure;
% subplot(2,2,1)
% %cspy(ConnMatrixPre(matIndex_commdet,matIndex_commdet),'Colormap',colorcet('R3'),'Levels',255,'MarkerSize',12);
% cspy(tempMat,'Colormap',colorcet('R3'),'Levels',255,'MarkerSize',12);
% 
% for i = 1:length(blocks)
%     line([0,size(tempMat,2)]+0.5,[blocks(i),blocks(i)]+0.5,'color','k');
%     line([blocks(i),blocks(i)]+0.5,[0,size(tempMat,2)]+0.5,'color','k');
% end
% 
% %  line([0,size(tempMat,2)]+0.5,[size(Vsort,1),size(Vsort,1)]+0.5,'color','k');
% %  line([size(Vsort,1),size(Vsort,1)]+0.5,[0,size(tempMat,2)]+0.5,'color','k');
% %  
% %  line([0,size(tempMat,2)]+0.5,[size(Vsort,1)+gap,size(Vsort,1)+gap]+0.5,'color','k');
% %  line([size(Vsort,1)+gap,size(Vsort,1)+gap]+0.5,[0,size(tempMat,2)]+0.5,'color','k');
% % % 
% %  line([0,size(tempMat,2)+0.5],[size(Vsort,1)+32+gap,size(Vsort,1)+32+gap]+0.5,'color','k');
% %  line([size(Vsort,1)+32+gap,size(Vsort,1)+32+gap]+0.5,[0,size(tempMat,2)+0.5],'color','k');
% 
% 
% set(gca,'XTick',[],'YTick',[]);
% box on;
% 
% locateRS = find(ismember(AllCellNames(idx),cellID(findRS)));
% text(locateRS,repmat(-20,1,numel(locateRS)), sprintf('\\downarrow'), 'HorizontalAlignment','center', 'FontWeight','bold','Color','k');


%% Pre and Post locations (long time)

% block1

block1_gamma038_WithNull.cellIDs(block1_gamma038_WithNull.cellIDs ==76202) = []; % remove mauthern for further analysis.
block1_gamma038_WithNull.dendSites = [];
block1_gamma038_WithNull.dendpartnerSites = [];
for i = 1:length(block1_gamma038_WithNull.cellIDs)
    [temp,preIDs] = SynapticPartners(block1_gamma038_WithNull.cellIDs(i),1,df);
    preIDs = preIDs(temp~=block1_gamma038_WithNull.cellIDs(i)); % remove self touches
    block1_gamma038_WithNull.dendSites = [block1_gamma038_WithNull.dendSites;PrePartnerCoordinates(preIDs,df)];
    block1_gamma038_WithNull.dendpartnerSites = [block1_gamma038_WithNull.dendpartnerSites;PostPartnerCoordinates(preIDs,df)];
    clear preIDs;
    clear temp;
    i
end
block1_gamma038_WithNull.dendSites = TransformPoints(block1_gamma038_WithNull.dendSites,0);
block1_gamma038_WithNull.dendpartnerSites = TransformPoints(block1_gamma038_WithNull.dendpartnerSites,0);


block2_gamma038_WithNull.dendSites = [];
block2_gamma038_WithNull.dendpartnerSites = [];
for i = 1:length(block2_gamma038_WithNull.cellIDs)
    [temp,preIDs] = SynapticPartners(block2_gamma038_WithNull.cellIDs(i),1,df);
    preIDs = preIDs(temp~=block2_gamma038_WithNull.cellIDs(i)); % remove self touches
    block2_gamma038_WithNull.dendSites = [block2_gamma038_WithNull.dendSites;PrePartnerCoordinates(preIDs,df)];
    block2_gamma038_WithNull.dendpartnerSites = [block2_gamma038_WithNull.dendpartnerSites;PostPartnerCoordinates(preIDs,df)];
    clear preIDs;
    i
end
block2_gamma038_WithNull.dendSites = TransformPoints(block2_gamma038_WithNull.dendSites,0);
block2_gamma038_WithNull.dendpartnerSites = TransformPoints(block2_gamma038_WithNull.dendpartnerSites,0);

% plot somata

block1_gamma038_WithNull.noSoma = [78158,80842,79560,78116,79503,79383,78125,77799,81147,77101,78121,78104,77239,81519,79771,79559,78914,79951,79240,80539,79054,81611,76936,79523,79407,79282,76626,77602,77607,81410,80472,79341,79732,77039,78119];
block2_gamma038_WithNull.noSoma = [77845,77152,79022,81295,79720,77151,78667,81395,77142,80821,81312,79062,80548,78633,80743,77636,80629,77461,77163,81336,80510,77797,77822,77656,76618,77389,77872,76697,76625,77434,77467,78351,78357,77621,76627,77868,77460,80679,77844,78650,77132,77126,81431,78406,78544,81423,77433,81637,80746,80625,81683,77805,79743,77821,79067,81792,76622,78404,78629,77465,78421,80757,79042,77162,79074,77122,79033,77651,77447,78255,80681,78543,78150];

block1_gamma038_WithNull.origins = getOrigin(block1_gamma038_WithNull.cellIDs);
block2_gamma038_WithNull.origins = getOrigin(block2_gamma038_WithNull.cellIDs);

block1_gamma038_WithNull.color = [1,0.5,0];
block2_gamma038_WithNull.color = [0,0.5,1];

figure;

% scatter plot for neurons with somata inside the volume.

 h = scatterhist([block1_gamma038_WithNull.origins(~ismember(block1_gamma038_WithNull.cellIDs,block1_gamma038_WithNull.noSoma),2);...
     block2_gamma038_WithNull.origins(~ismember(block2_gamma038_WithNull.cellIDs,block2_gamma038_WithNull.noSoma),2)],...
     [block1_gamma038_WithNull.origins(~ismember(block1_gamma038_WithNull.cellIDs,block1_gamma038_WithNull.noSoma),3);...
     block2_gamma038_WithNull.origins(~ismember(block2_gamma038_WithNull.cellIDs,block2_gamma038_WithNull.noSoma),3)],...
     'Group',[ones(sum(~ismember(block1_gamma038_WithNull.cellIDs,block1_gamma038_WithNull.noSoma)),1);2*ones(sum(~ismember(block2_gamma038_WithNull.cellIDs,block2_gamma038_WithNull.noSoma)),1)],...
     'Kernel','on','color',[block1_gamma038_WithNull.color;block2_gamma038_WithNull.color],'Marker','o','MarkerSize',5,'LineWidth',2);
  set(h(1),'YDir','reverse');
 set(h(3),'XDir','reverse');
daspect([1,1,1]);

% all somata scatter plot
figure

 h = scatterhist([block1_gamma038_WithNull.origins(:,1);block2_gamma038_WithNull.origins(:,1)],[block1_gamma038_WithNull.origins(:,2);block2_gamma038_WithNull.origins(:,2)],...
     'Group',[ones(length(block1_gamma038_WithNull.cellIDs),1);2*ones(length(block2_gamma038_WithNull.cellIDs),1)],...
     'Kernel','on','color',[block1_gamma038_WithNull.color;block2_gamma038_WithNull.color],'Marker','o','MarkerSize',5,'LineWidth',2);
set(h(1),'YDir','reverse','YLim',[400,800],'XLim',[250,400]);
set(h(2),'Visible','on');
set(h(3),'XDir','reverse','Visible','on');
daspect([1,1,1]);

% plot histograms of the inputs
Xcords_Pre = [block1_gamma038_WithNull.dendSites(1:10:end,1);block2_gamma038_WithNull.dendSites(1:10:end,1)];
Ycords_Pre = [block1_gamma038_WithNull.dendSites(1:10:end,2);block2_gamma038_WithNull.dendSites(1:10:end,2)];
Zcords_pre = [block1_gamma038_WithNull.dendSites(1:10:end,3);block2_gamma038_WithNull.dendSites(1:10:end,3)];

GpIDs_Pre = [ones(length(block1_gamma038_WithNull.dendSites(1:10:end,1)),1);2*ones(length(block2_gamma038_WithNull.dendSites(1:10:end,1)),1)];

h = scatterhist(Xcords_Pre,Ycords_Pre,'Group',GpIDs_Pre,'Kernel','on','color',[block1_gamma038_WithNull.color;block2_gamma038_WithNull.color],...
    'Marker','.','MarkerSize',1,'LineWidth',2);
set(h(1),'YDir','reverse','YLim',[400,800],'XLim',[250,400]);
set(h(2),'Visible','on');
set(h(3),'XDir','reverse','Visible','on');
daspect([1,1,1]);

figure;
h = scatterhist(Ycords_Pre,Zcords_pre,'Group',GpIDs_Pre,'Kernel','on','color',[block1_gamma038_WithNull.color;block2_gamma038_WithNull.color],...
    'Marker','.','MarkerSize',1,'LineWidth',2);
set(h(1),'YDir','reverse');
set(h(3),'XDir','reverse');
daspect([1,1,1]);




% n =100;
% Xlin = linspace(min(block1_gamma038_WithNull.dendSites(:,1)), max(block1_gamma038_WithNull.dendSites(:,1)),n);
% Ylin = linspace(min(block1_gamma038_WithNull.dendSites(:,2)), max(block1_gamma038_WithNull.dendSites(:,2)),n);
% 
% Xr = interp1(Xlin,1:numel(Xlin),block1_gamma038_WithNull.dendSites(:,1),'nearest');
% Yr = interp1(Ylin,1:numel(Ylin),block1_gamma038_WithNull.dendSites(:,2),'nearest');
% 
% Z = accumarray([Xr,Yr],1,[n,n]);
% figure;
% contour(Z,'color',col1);


% h = scatterhist(Ycords_Pre,Zcords_pre,'Group',GpIDs_Pre,'Kernel','on','color',[col1;col2],...
%     'Marker','.','MarkerSize',1,'LineWidth',2);

% block2_gamma038_WithNull
block1_gamma038_WithNull.axonSites = [];
block1_gamma038_WithNull.axonPartnerSites = [];
for i = 1:length(block1_gamma038_WithNull.cellIDs)
    [temp,preIDs] = SynapticPartners(block1_gamma038_WithNull.cellIDs(i),2,df);
    preIDs = preIDs(temp~=block1_gamma038_WithNull.cellIDs(i)); % remove self touches
    block1_gamma038_WithNull.axonSites = [block1_gamma038_WithNull.axonSites;PostPartnerCoordinates(preIDs,df)];
    block1_gamma038_WithNull.axonPartnerSites = [block1_gamma038_WithNull.axonPartnerSites;PrePartnerCoordinates(preIDs,df)];
    clear preIDs;
end
block1_gamma038_WithNull.axonSites = TransformPoints(block1_gamma038_WithNull.axonSites,0);
block1_gamma038_WithNull.axonPartnerSites = TransformPoints(block1_gamma038_WithNull.axonPartnerSites,0);

block2_gamma038_WithNull.axonSites = [];
block2_gamma038_WithNull.axonPartnerSites = [];
for i = 1:length(block2_gamma038_WithNull.cellIDs)
    [temp,preIDs] = SynapticPartners(block2_gamma038_WithNull.cellIDs(i),2,df);
    preIDs = preIDs(temp~=block2_gamma038_WithNull.cellIDs(i)); % remove self touches
    block2_gamma038_WithNull.axonSites = [block2_gamma038_WithNull.axonSites;PostPartnerCoordinates(preIDs,df)];
    block2_gamma038_WithNull.axonPartnerSites = [block2_gamma038_WithNull.axonPartnerSites;PrePartnerCoordinates(preIDs,df)];
    clear preIDs;
end
block2_gamma038_WithNull.axonSites  = TransformPoints(block2_gamma038_WithNull.axonSites ,0);
block2_gamma038_WithNull.axonPartnerSites = TransformPoints(block2_gamma038_WithNull.axonPartnerSites,0);

% plot histogram of the outputs
Xcords_Post = [block1_gamma038_WithNull.axonSites(1:10:end,1);block2_gamma038_WithNull.axonSites(1:10:end,1)];
Ycords_Post = [block1_gamma038_WithNull.axonSites(1:10:end,2);block2_gamma038_WithNull.axonSites(1:10:end,2)];
Zcords_Post = [block1_gamma038_WithNull.axonSites(1:10:end,3);block2_gamma038_WithNull.axonSites(1:10:end,3)];

figure;
 GpIDs_Post = [ones(length(block1_gamma038_WithNull.axonSites(1:10:end,1)),1);2*ones(length(block2_gamma038_WithNull.axonSites(1:10:end,1)),1)];
h = scatterhist(Xcords_Post,Ycords_Post,'Group',GpIDs_Post,'Kernel','on','color',[block1_gamma038_WithNull.color;block2_gamma038_WithNull.color],...
    'Marker','.','MarkerSize',0.5,'LineWidth',2);
set(h(1),'YDir','reverse','YLim',[400,800],'XLim',[250,400]);
set(h(2),'Visible','on');
set(h(3),'XDir','reverse','Visible','on');
daspect([1,1,1]);

figure;

h = scatterhist(Ycords_Post,Zcords_Post,'Group',GpIDs_Post,'Kernel','on','color',[block1_gamma038_WithNull.color;block2_gamma038_WithNull.color],...
    'Marker','.','MarkerSize',0.5,'LineWidth',2);
 set(h(1),'YDir','reverse');
 set(h(3),'XDir','reverse');
daspect([1,1,1]);


%% Peters Rule
% Do block1_gamma038_WithNull neurons have the potential to make synapses onto the other
% block.

% block1_gamma038_WithNull --> block1_gamma038_WithNull (actual synapses/ potential synapses)
% block1_gamma038_WithNull --> block2_gamma038_WithNull

for i = 1:length(block1_gamma038_WithNull.cellIDs)
    [prePartner,prePartnerPSD] = SynapticPartners(block1_gamma038_WithNull.cellIDs(i),1,df);
    % remove self touches
    block1_gamma038_WithNull.block1_gamma038_WithNulltoblock1_gamma038_WithNullActual(i) = sum(ismember(prePartner,setdiff(block1_gamma038_WithNull.cellIDs,block1_gamma038_WithNull.cellIDs(i)))); 
    prePartnerPSD = prePartnerPSD(prePartner~=block1_gamma038_WithNull.cellIDs(i));

    inputs = PrePartnerCoordinates(prePartnerPSD,df); % location of presynaptic partner, on the axon
    inputs = TransformPoints(inputs,0);
    % distance matrix block1_gamma038_WithNull -->block1_gamma038_WithNull
    distMat11 = pdist2(block1_gamma038_WithNull.axonPartnerSites,inputs);
    block1_gamma038_WithNull.actualSynapses11(i) = length(find(distMat11 == 0));
    block1_gamma038_WithNull.potentialSynpses11_05(i) = length(find(distMat11 <= 0.5));
    block1_gamma038_WithNull.potentialSynpses11_1(i) = length(find(distMat11 <= 1));
    block1_gamma038_WithNull.potentialSynpses11_2(i) = length(find(distMat11 <= 2));
    block1_gamma038_WithNull.potentialSynpses11_5(i) = length(find(distMat11 <= 5));
    block1_gamma038_WithNull.potentialSynpses11_10(i) = length(find(distMat11 <=10));


    %distance Matrix block1_gamma038_WithNull -->block2_gamma038_WithNull
    distMat12 = pdist2(block2_gamma038_WithNull.axonPartnerSites,inputs);
    block1_gamma038_WithNull.actualSynapses12(i) = length(find(distMat12 == 0));
    block1_gamma038_WithNull.potentialSynpses12_05(i) = length(find(distMat12 <= 0.5));
    block1_gamma038_WithNull.potentialSynpses12_1(i) = length(find(distMat12 <= 1));
    block1_gamma038_WithNull.potentialSynpses12_2(i) = length(find(distMat12 <= 2));
    block1_gamma038_WithNull.potentialSynpses12_5(i) = length(find(distMat12 <= 5));
    block1_gamma038_WithNull.potentialSynpses12_10(i) = length(find(distMat12 <= 10));

    
    clear PrePartner;
    clear inputs;
    clear distMat11;
    clear distMat12;
end

% block2_gamma038_WithNull --> block2_gamma038_WithNull (actual synapses/ potential synapses)
% block2_gamma038_WithNull --> block1_gamma038_WithNull

for i = 1:length(block2_gamma038_WithNull.cellIDs)
    [prePartner,prePartnerPSD] = SynapticPartners(block2_gamma038_WithNull.cellIDs(i),1,df);
    % remove self touches
    block2_gamma038_WithNull.block2_gamma038_WithNulltoblock2_gamma038_WithNullActual(i) = sum(ismember(prePartner,setdiff(block2_gamma038_WithNull.cellIDs,block2_gamma038_WithNull.cellIDs(i)))); 
    prePartnerPSD = prePartnerPSD(prePartner~=block2_gamma038_WithNull.cellIDs(i));
    inputs = PrePartnerCoordinates(prePartnerPSD,df); % location of presynaptic partner, on the axon
    inputs = TransformPoints(inputs,0);
    % distance matrix block2_gamma038_WithNull -->block2_gamma038_WithNull
    distMat22 = pdist2(block2_gamma038_WithNull.axonPartnerSites,inputs);
    block2_gamma038_WithNull.actualSynapses22(i) = length(find(distMat22 == 0));
    block2_gamma038_WithNull.potentialSynpses22_05(i) = length(find(distMat22 <= 0.5));
    block2_gamma038_WithNull.potentialSynpses22_1(i) = length(find(distMat22 <= 1));
    block2_gamma038_WithNull.potentialSynpses22_2(i) = length(find(distMat22 <= 2));
    block2_gamma038_WithNull.potentialSynpses22_5(i) = length(find(distMat22 <= 5));
    block2_gamma038_WithNull.potentialSynpses22_10(i) = length(find(distMat22 <= 10));


    %distance Matrix block2_gamma038_WithNull -->block1_gamma038_WithNull
    distMat21 = pdist2(block1_gamma038_WithNull.axonPartnerSites,inputs);
    block2_gamma038_WithNull.actualSynapses21(i) = length(find(distMat21 == 0));
    block2_gamma038_WithNull.potentialSynpses21_05(i) = length(find(distMat21 <= 0.5));
    block2_gamma038_WithNull.potentialSynpses21_1(i) = length(find(distMat21 <= 1));
    block2_gamma038_WithNull.potentialSynpses21_2(i) = length(find(distMat21 <= 2));
    block2_gamma038_WithNull.potentialSynpses21_5(i) = length(find(distMat21 <= 5));
    block2_gamma038_WithNull.potentialSynpses21_10(i) = length(find(distMat21 <= 10));

    
    clear PrePartner;
    clear inputs;
    clear distMat22;
    clear distMat21;
end

% figure;
% subplot(4,4,1)
% heatmap([mean(block1_gamma038_WithNull.actualSynapses11),mean( block1_gamma038_WithNull.actualSynapses12);mean( block2_gamma038_WithNull.actualSynapses21),mean(block2_gamma038_WithNull.actualSynapses22)],...
%     'ColorbarVisible','off');
% title('True synapses');
% 
% subplot(4,4,2)
% heatmap([mean(block1_gamma038_WithNull.potentialSynpses11_2),mean(block1_gamma038_WithNull.potentialSynpses12_2);mean( block2_gamma038_WithNull.potentialSynpses21_2),mean(block2_gamma038_WithNull.potentialSynpses22_2)],...
%     'ColorbarVisible','off');
% title('Potential synapses (2um)');
% 
% subplot(4,4,3)
% heatmap([mean(block1_gamma038_WithNull.potentialSynpses11_5),mean(block1_gamma038_WithNull.potentialSynpses12_5);mean( block2_gamma038_WithNull.potentialSynpses21_5),mean(block2_gamma038_WithNull.potentialSynpses22_5)],...
%     'ColorbarVisible','off');
% title('Potential synapses (5um)');
% 
% subplot(4,4,4)
% heatmap([mean(block1_gamma038_WithNull.potentialSynpses11_10),mean(block1_gamma038_WithNull.potentialSynpses12_10);mean( block2_gamma038_WithNull.potentialSynpses21_10),mean(block2_gamma038_WithNull.potentialSynpses22_10)],...
%     'ColorbarVisible','off');
% title('Potential synapses (10um)');


trueDiag = sum(block1_gamma038_WithNull.actualSynapses11)+sum(block2_gamma038_WithNull.actualSynapses22);
trueOffDiag = sum(block2_gamma038_WithNull.actualSynapses21)+sum(block1_gamma038_WithNull.actualSynapses12);

trueDiag_2 = sum(block1_gamma038_WithNull.potentialSynpses11_2)+sum(block2_gamma038_WithNull.potentialSynpses22_2);
trueOffDiag_2 = sum(block2_gamma038_WithNull.potentialSynpses21_2)+sum(block1_gamma038_WithNull.potentialSynpses12_2);

trueDiag_5 = sum(block1_gamma038_WithNull.potentialSynpses11_5)+sum(block2_gamma038_WithNull.potentialSynpses22_5);
trueOffDiag_5 = sum(block2_gamma038_WithNull.potentialSynpses21_5)+sum(block1_gamma038_WithNull.potentialSynpses12_5);

trueDiag_10 = sum(block1_gamma038_WithNull.potentialSynpses11_10)+sum(block2_gamma038_WithNull.potentialSynpses22_10);
trueOffDiag_10 = sum(block2_gamma038_WithNull.potentialSynpses21_10)+sum(block1_gamma038_WithNull.potentialSynpses12_10);

figure;
subplot(4,4,1)
plot([1,2,3,4],[trueDiag/trueOffDiag,trueDiag_2/trueOffDiag_2,trueDiag_5/trueOffDiag_5,trueDiag_10/trueOffDiag_10],'-ko','LineWidth',2);
set(gca,'XTick',[1:4],'XTickLabels',[{'true'},2,5,10],'YLim',[1,6]);
box off;
daspect([1,1,1]);
offsetAxes(gca);
xlabel('radius (\mum)');
ylabel('diagonal sum/off-diagonal sum')

% subplot(4,4,5)
% % scatter(block1_gamma038_WithNull.actualSynapses11,block1_gamma038_WithNull.potentialSynpses11_5);
% hold on
% % scatter(block1_gamma038_WithNull.actualSynapses12,block1_gamma038_WithNull.potentialSynpses12_5);
% % scatter(block2_gamma038_WithNull.actualSynapses22,block2_gamma038_WithNull.potentialSynpses22_5);
% % scatter(block2_gamma038_WithNull.actualSynapses21,block2_gamma038_WithNull.potentialSynpses21_5);
% %line([0,100],[0,100],'color','k');
% showfit(ezfit(block1_gamma038_WithNull.actualSynapses11,block1_gamma038_WithNull.potentialSynpses11_5,'affine'),'fitcolor',col1,'dispeqboxmode','off');
% showfit(ezfit(block1_gamma038_WithNull.actualSynapses12,block1_gamma038_WithNull.potentialSynpses12_5,'affine'),'fitcolor',col1,'dispeqboxmode','off','fitlinestyle',':');
% showfit(ezfit(block2_gamma038_WithNull.actualSynapses22,block2_gamma038_WithNull.potentialSynpses22_5,'affine'),'fitcolor',col2,'dispeqboxmode','off');
% showfit(ezfit(block2_gamma038_WithNull.actualSynapses21,block2_gamma038_WithNull.potentialSynpses21_5,'affine'),'fitcolor',col2,'dispeqboxmode','off','fitlinestyle',':');
% set(gca,'XLim',[0,100],'YLim',[0,100000]);
% legend({'b1-->b1','b1-->b2','b2-->b2','b2-->b1'});
% %axis square;
% xlabel('True Synapses');
% ylabel('Potential Synapses');

% peters analysis for block1_gamma038_WithNull-->ABD and block2_gamma038_WithNull -->ABD

ABDr_CellIDs = [82140,82145,82143,77648,82146,77302,77301,82194,77705,77710,77305,77672,77300,77709,77661];
ABDc_CellIDs = [82213,77654,77646,82212,77295,81172,77657,82197,77682,77154, 77652,77658,77628,77292,77688,82195,77296];
ABDIr_CellIDs = [78553,77150,77886,78547,77631,77158,77665,78552,77668,77618,77634,78556];
ABDIc_CellIDs = [78574,79051,77148,77625,77144,77692,77641,79066,77640,77643];

% get all axonal sites for ABDm neruons

motorNeurons = [ABDr_CellIDs,ABDc_CellIDs];
motorNeurons_dendPartners = [];
for i = 1:length(motorNeurons)
    [temp,prePSDid] = SynapticPartners(motorNeurons(i),1,df);
    prePSDid = prePSDid(temp~=prePSDid(i)); % remove self touches
    motorNeurons_dendPartners = [motorNeurons_dendPartners;PrePartnerCoordinates(prePSDid,df)];
    clear temp;
    clear prePSDid;
end
motorNeurons_dendPartners = TransformPoints(motorNeurons_dendPartners,0);

% interneuron dendrite partner sites
interNeurons = [ABDIr_CellIDs,ABDIc_CellIDs];
interNeurons_dendPartners = [];
for i = 1:length(interNeurons)
    [temp,prePSDid] = SynapticPartners(interNeurons(i),1,df);
    prePSDid = prePSDid(temp~=prePSDid(i)); % remove self touches
    interNeurons_dendPartners = [interNeurons_dendPartners;PrePartnerCoordinates(prePSDid,df)];
    clear temp;
    clear prePSDid;
end
interNeurons_dendPartners = TransformPoints(interNeurons_dendPartners,0);



for i = 1:length(block1_gamma038_WithNull.cellIDs)
    [temp,postID] = SynapticPartners(block1_gamma038_WithNull.cellIDs(i),2,df);
    postID = postID(temp~=block1_gamma038_WithNull.cellIDs(i));
    temp2= PrePartnerCoordinates(postID,df);
    temp2 = TransformPoints(temp2,0);
    
    distMatmotor = pdist2(temp2,motorNeurons_dendPartners);
    distMatinter = pdist2(temp2,interNeurons_dendPartners);
    
    block1_gamma038_WithNull.actualSynapsesMotor(i) = length(find(distMatmotor == 0));
    block1_gamma038_WithNull.potentialSynapsesMotor_2(i) = length(find(distMatmotor <=2));
    block1_gamma038_WithNull.potentialSynapsesMotor_5(i) = length(find(distMatmotor <=5));
    block1_gamma038_WithNull.potentialSynapsesMotor_10(i) = length(find(distMatmotor <=10));
    
    block1_gamma038_WithNull.actualSynapsesInter(i) = length(find(distMatinter == 0));
    block1_gamma038_WithNull.potentialSynapsesInter_2(i) = length(find(distMatinter <=2));
    block1_gamma038_WithNull.potentialSynapsesInter_5(i) = length(find(distMatinter <=5));
    block1_gamma038_WithNull.potentialSynapsesInter_10(i) = length(find(distMatinter <=10));
    
    %block1_gamma038_WithNull.axonalSites = [block1_gamma038_WithNull.axonalSites;temp2]
    
    clear temp
    clear postID;
    clear distMatmotor
    clear distMatinter
end


for i = 1:length(block2_gamma038_WithNull.cellIDs)
    [temp,postID] = SynapticPartners(block2_gamma038_WithNull.cellIDs(i),2,df);
    postID = postID(temp~=block2_gamma038_WithNull.cellIDs(i));
    temp2= PrePartnerCoordinates(postID,df);
    temp2 = TransformPoints(temp2,0);
    
    distMatmotor = pdist2(temp2,motorNeurons_dendPartners);
    distMatinter = pdist2(temp2,interNeurons_dendPartners);
    
    block2_gamma038_WithNull.actualSynapsesMotor(i) = length(find(distMatmotor == 0));
    block2_gamma038_WithNull.potentialSynapsesMotor_2(i) = length(find(distMatmotor <=2));
    block2_gamma038_WithNull.potentialSynapsesMotor_5(i) = length(find(distMatmotor <=5));
    block2_gamma038_WithNull.potentialSynapsesMotor_10(i) = length(find(distMatmotor <=10));
    
    block2_gamma038_WithNull.actualSynapsesInter(i) = length(find(distMatinter == 0));
    block2_gamma038_WithNull.potentialSynapsesInter_2(i) = length(find(distMatinter <=2));
    block2_gamma038_WithNull.potentialSynapsesInter_5(i) = length(find(distMatinter <=5));
    block2_gamma038_WithNull.potentialSynapsesInter_10(i) = length(find(distMatinter <=10));
    
    %block2_gamma038_WithNull.axonalSites = [block2_gamma038_WithNull.axonalSites;temp2]
    
    clear temp
    clear postID;
    clear distMatmotor
    clear distMatinter
end


% Potential synapse, mod to ABD

mod1ABD = sum(block1_gamma038_WithNull.actualSynapsesMotor)+sum(block1_gamma038_WithNull.actualSynapsesInter);
mod2ABD = sum(block2_gamma038_WithNull.actualSynapsesMotor)+sum(block2_gamma038_WithNull.actualSynapsesInter);

mod1ABD_2 = sum(block1_gamma038_WithNull.potentialSynapsesMotor_2)+sum(block1_gamma038_WithNull.potentialSynapsesInter_2);
mod1ABD_5 = sum(block1_gamma038_WithNull.potentialSynapsesMotor_5)+sum(block1_gamma038_WithNull.potentialSynapsesInter_5);
mod1ABD_10 = sum(block1_gamma038_WithNull.potentialSynapsesMotor_10)+sum(block1_gamma038_WithNull.potentialSynapsesInter_10);

mod2ABD_2 = sum(block2_gamma038_WithNull.potentialSynapsesMotor_2)+sum(block2_gamma038_WithNull.potentialSynapsesInter_2);
mod2ABD_5 = sum(block2_gamma038_WithNull.potentialSynapsesMotor_5)+sum(block2_gamma038_WithNull.potentialSynapsesInter_5);
mod2ABD_10 = sum(block2_gamma038_WithNull.potentialSynapsesMotor_10)+sum(block2_gamma038_WithNull.potentialSynapsesInter_10);

subplot(4,4,6)
% semilogy([mod1ABD,mod1ABD_2,mod1ABD_5,mod1ABD_10],'-o','color',[1,0.5,0],'LineWidth',2);
% hold on
% semilogy([mod2ABD,mod2ABD_2,mod2ABD_5,mod2ABD_10],'-o','color',[0,0.5,1],'LineWidth',2);

% plot([1,2,3,4],[mod1ABD,mod1ABD_2,mod1ABD_5,mod1ABD_10]./[mod2ABD,mod2ABD_2,mod2ABD_5,mod2ABD_10],...
%     '-ok','LineWidth',2);
% yyaxis right
plot([1,2,3,4],[mod2ABD,mod2ABD_2,mod2ABD_5,mod2ABD_10]./[mod1ABD,mod1ABD_2,mod1ABD_5,mod1ABD_10],...
    '-ok','LineWidth',2);
%daspect([1,2,1]);
box off;
set(gca,'XTick',[1:4],'XTickLabels',[{'true'},2,5,10],'YLim',[1,16],'YTick',[1:3:16]);
ylabel('mod1/mod2 synapses');
xlabel('radius (\mum)');
daspect([1,3,1]);
offsetAxes(gca);

%axis square;
%legend({'mod1','mod2'},'Location','bestoutside');

%% calculate recurrent fractions

for i = 1:length(block2_gamma038_WithNull.cellIDs)
    A = SynapticPartners(block2_gamma038_WithNull.cellIDs(i),1,df);
    block2_gamma038_WithNull.recurrentFraction(i) = sum(ismember(A,block2_gamma038_WithNull.cellIDs))./length(A);
    clear A;
end

for i = 1:length(block1_gamma038_WithNull.cellIDs)
    A = SynapticPartners(block1_gamma038_WithNull.cellIDs(i),1,df);
    block1_gamma038_WithNull.recurrentFraction(i) = sum(ismember(A,block1_gamma038_WithNull.cellIDs))./length(A);
    clear A;
end

% submodules
for i = 1:length(block1_subMod_gamm075.cellID)
    A = SynapticPartners(block1_subMod_gamm075.cellID(i),1,df);
    block1_subMod_gamm075.recurrentFraction(i) = sum(ismember(A,block1_subMod_gamm075.cellID))./length(A);
    clear A;
end

for i = 1:length(block2_subMod_gamm075.cellID)
    A = SynapticPartners(block2_subMod_gamm075.cellID(i),1,df);
    block2_subMod_gamm075.recurrentFraction(i) = sum(ismember(A,block2_subMod_gamm075.cellID))./length(A);
    clear A;
end


% feedforward 
%%
figure;
subplot(4,4,1)
histogram(rmoutliers(block2_gamma038_WithNull.recurrentFraction),'BinWidth',0.05,'EdgeColor',[0,0.5,1],'DisplayStyle','stairs','LineWidth',2);
hold on;
histogram(rmoutliers(block1_gamma038_WithNull.recurrentFraction),'BinWidth',0.05,'EdgeColor',[1,0.5,0],'DisplayStyle','stairs','LineWidth',2);
legend({'mod2','mod1'},'Location','bestoutside');
axis square;
xlabel(' Recurrent fraction');
ylabel(' count');
box off;
offsetAxes(gca);

block1_subMod_gamm075.color = [0.6207    0.3103    0.2759];
block2_subMod_gamm075.color = [0.5172    0.5172    1.0000];


subplot(4,4,2)
histogram(rmoutliers(block1_subMod_gamm075.recurrentFraction),'BinWidth',0.05,'EdgeColor',block1_subMod_gamm075.color,'DisplayStyle','stairs','LineWidth',2);
hold on;
histogram(rmoutliers(block2_subMod_gamm075.recurrentFraction),'BinWidth',0.05,'EdgeColor',block2_subMod_gamm075.color,'DisplayStyle','stairs','LineWidth',2);
legend({'mod2a','mod2b'},'Location','bestoutside');
axis square;
xlabel(' Recurrent fraction');
ylabel(' count');
box off;
offsetAxes(gca);

%%
%cbl = [80478,77251,80529,80341,78236,80511,80482,80502,80552,80514,80554,80352,80506,80496,80519,80490,80518,80339,78239,80524,78280,80547];

% 
% matOrder_commdet = [DOs,...
%     block1_gamma038_WithNull,block2_gamma038_WithNull,Block3,...
%     ABDr_CellIDs,ABDc_CellIDs,ABDIr_CellIDs,ABDIc_CellIDs];
%    
% [~,matIndex_commdet] = ismember(matOrder_commdet,AllCells);
% 
% connMatForPlot = ConnMatrixPre(matIndex_commdet,matIndex_commdet);
% 
% for i = 1:size(matOrder_commdet,2)
%     matOrder_commdet_totalInputs(i) = length(SynapticPartners(matOrder_commdet(i),1,df));
%     connMatForPlot_Norm(i,:) = connMatForPlot(i,:)./matOrder_commdet_totalInputs(i);
% end
% connMatForPlot_Norm(connMatForPlot_Norm>0.05) = 0.05;

%% plot loation of neurons in block5 and block6

% transform_swc_AV(block5_gamma09.cellIDs,block5_gamma09.color,[],true,false);
% transform_swc_AV(block6_gamma09.cellIDs,block6_gamma09.color,[],true,false);

transform_swc_AV(block1_subMod_gamm075.cellID,hex2rgb('#cc5643'),[],true,false);
transform_swc_AV(block2_subMod_gamm075.cellID,hex2rgb('#8773c9'),[],false,false);

% get origins

block1_subMod_gamm075.Origins = getOrigin(block1_subMod_gamm075.cellID);
block2_subMod_gamm075.Origins = getOrigin(block2_subMod_gamm075.cellID);

% scatter histogram of spatial locations

Xcords_Pre = [block1_subMod_gamm075.Origins(:,1);block2_subMod_gamm075.Origins(:,1)];
Ycords_Pre = [block1_subMod_gamm075.Origins(:,2);block2_subMod_gamm075.Origins(:,2)];
Zcords_pre = [block1_subMod_gamm075.Origins(:,3);block2_subMod_gamm075.Origins(:,3)];

GpIDs_Pre = [ones(length(block1_subMod_gamm075.Origins(:,1)),1); 2*ones(length(block2_subMod_gamm075.Origins(:,1)),1)];

h = scatterhist(Xcords_Pre,Ycords_Pre,'Group',GpIDs_Pre,'Kernel','on','color',[block1_subMod_gamm075.color;block2_subMod_gamm075.color],...
    'Marker','o','MarkerSize',5,'LineWidth',2);
set(h(1),'YDir','reverse','YLim',[400,800],'XLim',[250,400]);
set(h(3),'XDir','reverse');
daspect([1,1,1]);

% location of presynapsess

block1_subMod_gamm075.dendSites = [];
block1_subMod_gamm075.dendpartnerSites = [];
for i = 1:length(block1_subMod_gamm075.cellID)
    [temp,preIDs] = SynapticPartners(block1_subMod_gamm075.cellID(i),1,df);
    preIDs = preIDs(temp~=block1_subMod_gamm075.cellID(i)); % remove self touches
    block1_subMod_gamm075.dendSites = [block1_subMod_gamm075.dendSites;PrePartnerCoordinates(preIDs,df)];
    block1_subMod_gamm075.dendpartnerSites = [block1_subMod_gamm075.dendpartnerSites;PostPartnerCoordinates(preIDs,df)];
    clear preIDs;
    clear temp;
end
block1_subMod_gamm075.dendSites = TransformPoints(block1_subMod_gamm075.dendSites,0);
block1_subMod_gamm075.dendpartnerSites = TransformPoints(block1_subMod_gamm075.dendpartnerSites,0);


block2_subMod_gamm075.dendSites = [];
block2_subMod_gamm075.dendpartnerSites = [];
for i = 1:length(block2_subMod_gamm075.cellID)
    [temp,preIDs] = SynapticPartners(block2_subMod_gamm075.cellID(i),1,df);
    preIDs = preIDs(temp~=block2_subMod_gamm075.cellID(i)); % remove self touches
    block2_subMod_gamm075.dendSites = [block2_subMod_gamm075.dendSites;PrePartnerCoordinates(preIDs,df)];
    block2_subMod_gamm075.dendpartnerSites = [block2_subMod_gamm075.dendpartnerSites;PostPartnerCoordinates(preIDs,df)];
    clear preIDs;
    clear temp;
end
block2_subMod_gamm075.dendSites = TransformPoints(block2_subMod_gamm075.dendSites,0);
block2_subMod_gamm075.dendpartnerSites = TransformPoints(block2_subMod_gamm075.dendpartnerSites,0);


% location of postsynapses
block1_subMod_gamm075.axonSites = [];
block1_subMod_gamm075.axonPartnerSites = [];
for i = 1:length(block1_subMod_gamm075.cellID)
    [temp,preIDs] = SynapticPartners(block1_subMod_gamm075.cellID(i),2,df);
    preIDs = preIDs(temp~=block1_subMod_gamm075.cellID(i)); % remove self touches
    block1_subMod_gamm075.axonSites = [block1_subMod_gamm075.axonSites;PostPartnerCoordinates(preIDs,df)];
    block1_subMod_gamm075.axonPartnerSites = [block1_subMod_gamm075.axonPartnerSites;PrePartnerCoordinates(preIDs,df)];
    clear preIDs;
end
block1_subMod_gamm075.axonSites = TransformPoints(block1_subMod_gamm075.axonSites,0);
block1_subMod_gamm075.axonPartnerSites = TransformPoints(block1_subMod_gamm075.axonPartnerSites,0);

block2_subMod_gamm075.axonSites = [];
block2_subMod_gamm075.axonPartnerSites = [];
for i = 1:length(block2_subMod_gamm075.cellID)
    [temp,preIDs] = SynapticPartners(block2_subMod_gamm075.cellID(i),2,df);
    preIDs = preIDs(temp~=block2_subMod_gamm075.cellID(i)); % remove self touches
    block2_subMod_gamm075.axonSites = [block2_subMod_gamm075.axonSites;PostPartnerCoordinates(preIDs,df)];
    block2_subMod_gamm075.axonPartnerSites = [block2_subMod_gamm075.axonPartnerSites;PrePartnerCoordinates(preIDs,df)];
    clear preIDs;
end
block2_subMod_gamm075.axonSites = TransformPoints(block2_subMod_gamm075.axonSites,0);
block2_subMod_gamm075.axonPartnerSites = TransformPoints(block2_subMod_gamm075.axonPartnerSites,0);


% plot histograms of the inputs
Xcords_Pre = [block1_subMod_gamm075.dendSites(1:5:end,1);block2_subMod_gamm075.dendSites(1:5:end,1)];
Ycords_Pre = [block1_subMod_gamm075.dendSites(1:5:end,2);block2_subMod_gamm075.dendSites(1:5:end,2)];
Zcords_pre = [block1_subMod_gamm075.dendSites(1:5:end,3);block2_subMod_gamm075.dendSites(1:5:end,3)];

GpIDs_Pre = [ones(length(block1_subMod_gamm075.dendSites(1:5:end,1)),1);2*ones(length(block2_subMod_gamm075.dendSites(1:5:end,1)),1)];

h = scatterhist(Xcords_Pre,Ycords_Pre,'Group',GpIDs_Pre,'Kernel','on','color',[block1_subMod_gamm075.color;block2_subMod_gamm075.color],...
    'Marker','.','MarkerSize',1,'LineWidth',2);
set(h(1),'YDir','reverse','YLim',[400,800],'XLim',[250,400]);
set(h(3),'XDir','reverse');
daspect([1,1,1]);

% plot histogram of the outputs
Xcords_Post = [block1_subMod_gamm075.axonSites(1:5:end,1);block2_subMod_gamm075.axonSites(1:5:end,1)];
Ycords_Post = [block1_subMod_gamm075.axonSites(1:5:end,2);block2_subMod_gamm075.axonSites(1:5:end,2)];
Zcords_Post = [block1_subMod_gamm075.axonSites(1:5:end,3);block2_subMod_gamm075.axonSites(1:5:end,3)];

 GpIDs_Post = [ones(length(block1_subMod_gamm075.axonSites(1:5:end,1)),1);2*ones(length(block2_subMod_gamm075.axonSites(1:5:end,1)),1)];
h = scatterhist(Xcords_Post,Ycords_Post,'Group',GpIDs_Post,'Kernel','on','color',[block1_subMod_gamm075.color;block2_subMod_gamm075.color],...
    'Marker','.','MarkerSize',0.5,'LineWidth',2);
set(h(1),'YDir','reverse','YLim',[400,800],'XLim',[250,400]);
set(h(3),'XDir','reverse');
daspect([1,1,1]);

%% Peters analysis for submodules; block5 -->ABDm, block6 -->ABDi

% block1_ABDi --> block1_ABDi
for i = 1:length(block1_subMod_gamm075.cellID)
    [prePartner,prePartnerPSD] = SynapticPartners(block1_subMod_gamm075.cellID(i),1,df);
    % remove self touches
    %block1_gamma038_WithNull.block1_gamma038_WithNulltoblock1_gamma038_WithNullActual(i) = sum(ismember(prePartner,setdiff(block1_gamma038_WithNull.cellIDs,block1_gamma038_WithNull.cellIDs(i)))); 
    prePartnerPSD = prePartnerPSD(prePartner~=block1_subMod_gamm075.cellID(i));

    inputs = PrePartnerCoordinates(prePartnerPSD,df); % location of presynaptic partner, on the axon
    inputs = TransformPoints(inputs,0);
    % distance matrix block1_gamma038_WithNull -->block1_gamma038_WithNull
    distMat11 = pdist2(block1_subMod_gamm075.axonPartnerSites,inputs);
    block1_subMod_gamm075.actualSynapses11(i) = length(find(distMat11 == 0));
    block1_subMod_gamm075.potentialSynpses11_05(i) = length(find(distMat11 <= 0.5));
    block1_subMod_gamm075.potentialSynpses11_1(i) = length(find(distMat11 <= 1));
    block1_subMod_gamm075.potentialSynpses11_2(i) = length(find(distMat11 <= 2));
    block1_subMod_gamm075.potentialSynpses11_5(i) = length(find(distMat11 <= 5));
    block1_subMod_gamm075.potentialSynpses11_10(i) = length(find(distMat11 <=10));


    %distance Matrix block1_gamma038_WithNull -->block2_gamma038_WithNull
    distMat12 = pdist2(block2_subMod_gamm075.axonPartnerSites,inputs);
    block1_subMod_gamm075.actualSynapses12(i) = length(find(distMat12 == 0));
    block1_subMod_gamm075.potentialSynpses12_05(i) = length(find(distMat12 <= 0.5));
    block1_subMod_gamm075.potentialSynpses12_1(i) = length(find(distMat12 <= 1));
    block1_subMod_gamm075.potentialSynpses12_2(i) = length(find(distMat12 <= 2));
    block1_subMod_gamm075.potentialSynpses12_5(i) = length(find(distMat12 <= 5));
    block1_subMod_gamm075.potentialSynpses12_10(i) = length(find(distMat12 <= 10));

    
    clear PrePartner;
    clear inputs;
    clear distMat11;
    clear distMat12;
end

% block2_ABDm --> block2_ABDm
for i = 1:length(block2_subMod_gamm075.cellID)
    [prePartner,prePartnerPSD] = SynapticPartners(block2_subMod_gamm075.cellID(i),1,df);
    % remove self touches
    %block1_gamma038_WithNull.block1_gamma038_WithNulltoblock1_gamma038_WithNullActual(i) = sum(ismember(prePartner,setdiff(block1_gamma038_WithNull.cellIDs,block1_gamma038_WithNull.cellIDs(i)))); 
    prePartnerPSD = prePartnerPSD(prePartner~=block2_subMod_gamm075.cellID(i));

    inputs = PrePartnerCoordinates(prePartnerPSD,df); % location of presynaptic partner, on the axon
    inputs = TransformPoints(inputs,0);
    % distance matrix block1_gamma038_WithNull -->block1_gamma038_WithNull
    distMat22 = pdist2(block2_subMod_gamm075.axonPartnerSites,inputs);
    block2_subMod_gamm075.actualSynapses22(i) = length(find(distMat22 == 0));
    block2_subMod_gamm075.potentialSynpses22_05(i) = length(find(distMat22 <= 0.5));
    block2_subMod_gamm075.potentialSynpses22_1(i) = length(find(distMat22 <= 1));
    block2_subMod_gamm075.potentialSynpses22_2(i) = length(find(distMat22 <= 2));
    block2_subMod_gamm075.potentialSynpses22_5(i) = length(find(distMat22 <= 5));
    block2_subMod_gamm075.potentialSynpses22_10(i) = length(find(distMat22 <=10));


    %distance Matrix block1_gamma038_WithNull -->block2_gamma038_WithNull
    distMat21 = pdist2(block1_subMod_gamm075.axonPartnerSites,inputs);
    block2_subMod_gamm075.actualSynapses21(i) = length(find(distMat21 == 0));
    block2_subMod_gamm075.potentialSynpses21_05(i) = length(find(distMat21 <= 0.5));
    block2_subMod_gamm075.potentialSynpses21_1(i) = length(find(distMat21 <= 1));
    block2_subMod_gamm075.potentialSynpses21_2(i) = length(find(distMat21 <= 2));
    block2_subMod_gamm075.potentialSynpses21_5(i) = length(find(distMat21 <= 5));
    block2_subMod_gamm075.potentialSynpses21_10(i) = length(find(distMat21 <= 10));

    
    clear PrePartner;
    clear inputs;
    clear distMat22;
    clear distMat21;
end


trueDiag_submods = sum(block1_subMod_gamm075.actualSynapses11)+sum(block2_subMod_gamm075.actualSynapses22);
trueOffDiag_submods = sum(block1_subMod_gamm075.actualSynapses12) + sum( block2_subMod_gamm075.actualSynapses21);

trueDiag_submods_2 = sum(block1_subMod_gamm075.potentialSynpses11_2)+sum(block2_subMod_gamm075.potentialSynpses22_2);
trueOffDiag_submods_2 = sum(block1_subMod_gamm075.potentialSynpses12_2) + sum( block2_subMod_gamm075.potentialSynpses21_2);

trueDiag_submods_5 = sum(block1_subMod_gamm075.potentialSynpses11_5)+sum(block2_subMod_gamm075.potentialSynpses22_5);
trueOffDiag_submods_5 = sum(block1_subMod_gamm075.potentialSynpses12_5) + sum( block2_subMod_gamm075.potentialSynpses21_5);

trueDiag_submods_10 = sum(block1_subMod_gamm075.potentialSynpses11_10)+sum(block2_subMod_gamm075.potentialSynpses22_10);
trueOffDiag_submods_10 = sum(block1_subMod_gamm075.potentialSynpses12_10) + sum( block2_subMod_gamm075.potentialSynpses21_10);


figure;

subplot(4,4,1)

plot([1,2,3,4],[trueDiag_submods/trueOffDiag_submods,trueDiag_submods_2/trueOffDiag_submods_2,...
    trueDiag_submods_5/trueOffDiag_submods_5,trueDiag_submods_10/trueOffDiag_submods_10],'-ok','LineWidth',2);
set(gca,'XTick',[1:4],'XTickLabels',[{'true'},2,5,10],'YLim',[1,6]);
box off;
daspect([1,1,1]);
xlabel('radius (\mum)');
ylabel('$\displaystyle\frac{sum(Diag.)}{sum(off-Diag.)}$','Interpreter','latex')

% get all block1_ABDi-->motor

block1_subMod_gamm075.axonalSites = [];
for i = 1:length(block1_subMod_gamm075.cellID)
    [temp,postID] = SynapticPartners(block1_subMod_gamm075.cellID(i),2,df);
    postID = postID(temp~=block1_subMod_gamm075.cellID(i));
    temp2= PrePartnerCoordinates(postID,df);
    temp2 = TransformPoints(temp2,0);
    
    distMatmotor = pdist2(temp2,motorNeurons_dendPartners);
    distMatinter = pdist2(temp2,interNeurons_dendPartners);
    
    block1_subMod_gamm075.actualSynapsesMotor(i) = length(find(distMatmotor == 0));
    block1_subMod_gamm075.potentialSynapsesMotor_2(i) = length(find(distMatmotor <=2));
    block1_subMod_gamm075.potentialSynapsesMotor_5(i) = length(find(distMatmotor <=5));
    block1_subMod_gamm075.potentialSynapsesMotor_10(i) = length(find(distMatmotor <=10));
    
    block1_subMod_gamm075.actualSynapsesInter(i) = length(find(distMatinter == 0));
    block1_subMod_gamm075.potentialSynapsesInter_2(i) = length(find(distMatinter <=2));
    block1_subMod_gamm075.potentialSynapsesInter_5(i) = length(find(distMatinter <=5));
    block1_subMod_gamm075.potentialSynapsesInter_10(i) = length(find(distMatinter <=10));
    
    
    block1_subMod_gamm075.axonalSites = [block1_subMod_gamm075.axonalSites;temp2];
    
    clear temp
    clear postID;
end


% block6 --> motor

block2_subMod_gamm075.axonalSites = [];
for i = 1:length(block2_subMod_gamm075.cellID)
    [temp,postID] = SynapticPartners(block2_subMod_gamm075.cellID(i),2,df);
    postID = postID(temp~=block2_subMod_gamm075.cellID(i));
    temp2= PrePartnerCoordinates(postID,df);
    temp2 = TransformPoints(temp2,0);
    
    distMatmotor = pdist2(temp2,motorNeurons_dendPartners);
    distMatinter = pdist2(temp2,interNeurons_dendPartners);
    
    block2_subMod_gamm075.actualSynapsesMotor(i) = length(find(distMatmotor == 0));
    block2_subMod_gamm075.potentialSynapsesMotor_2(i) = length(find(distMatmotor <=2));
    block2_subMod_gamm075.potentialSynapsesMotor_5(i) = length(find(distMatmotor <=5));
    block2_subMod_gamm075.potentialSynapsesMotor_10(i) = length(find(distMatmotor <=10));
    
    block2_subMod_gamm075.actualSynapsesInter(i) = length(find(distMatinter == 0));
    block2_subMod_gamm075.potentialSynapsesInter_2(i) = length(find(distMatinter <=2));
    block2_subMod_gamm075.potentialSynapsesInter_5(i) = length(find(distMatinter <=5));
    block2_subMod_gamm075.potentialSynapsesInter_10(i) = length(find(distMatinter <=10));
    
    
    block2_subMod_gamm075.axonalSites = [block2_subMod_gamm075.axonalSites;temp2];
    
    
    clear temp
    clear postID;
end


subplot(4,4,2)
plot([1,2,3,4],[sum(block2_subMod_gamm075.actualSynapsesMotor),sum(block2_subMod_gamm075.potentialSynapsesMotor_2),...
    sum(block2_subMod_gamm075.potentialSynapsesMotor_5),sum(block2_subMod_gamm075.potentialSynapsesMotor_10)]./...
    [sum(block1_subMod_gamm075.actualSynapsesMotor),sum(block1_subMod_gamm075.potentialSynapsesMotor_2),...
    sum(block1_subMod_gamm075.potentialSynapsesMotor_5),sum(block1_subMod_gamm075.potentialSynapsesMotor_10)],...
    '-ok','LineWidth',2);
set(gca,'XTick',[1:4],'XTickLabels',[{'true'},2,5,10]);
%ylabel('mod2b/mod2a on ABD_M');
box off;
offsetAxes(gca);
hold on;
%subplot(4,4,3)
plot([1,2,3,4],[sum(block1_subMod_gamm075.actualSynapsesInter),sum(block1_subMod_gamm075.potentialSynapsesInter_2),...
    sum(block1_subMod_gamm075.potentialSynapsesInter_5),sum(block1_subMod_gamm075.potentialSynapsesInter_10)]./...
    [sum(block2_subMod_gamm075.actualSynapsesInter),sum(block2_subMod_gamm075.potentialSynapsesInter_2),...
    sum(block2_subMod_gamm075.potentialSynapsesInter_5),sum(block2_subMod_gamm075.potentialSynapsesInter_10)],...
    '-or','LineWidth',2);
%set(gca,'XTick',[1:4],'XTickLabels',[{'true'},2,5,10]);
%ylabel('mod2a/mod2b onto ABD_I');
legend({'mod2b/mod2a -> ABD_M','mod2a/mod2b -> ABD_I'})
box off;
%offsetAxes(gca);

subplot(4,4,3)
plot([1,2,3,4],[(sum(block1_subMod_gamm075.actualSynapsesInter)+sum(block2_subMod_gamm075.actualSynapsesMotor))/(sum(block1_subMod_gamm075.actualSynapsesMotor)+sum( block2_subMod_gamm075.actualSynapsesInter));...
                (sum(block2_subMod_gamm075.potentialSynapsesInter_2)+sum(block2_subMod_gamm075.potentialSynapsesMotor_2))/(sum(block1_subMod_gamm075.potentialSynapsesMotor_2)+sum(block2_subMod_gamm075.potentialSynapsesInter_2));...
                (sum(block2_subMod_gamm075.potentialSynapsesInter_5)+sum(block2_subMod_gamm075.potentialSynapsesMotor_5))/(sum(block1_subMod_gamm075.potentialSynapsesMotor_5)+sum(block2_subMod_gamm075.potentialSynapsesInter_5));...
                (sum(block2_subMod_gamm075.potentialSynapsesInter_10)+sum(block2_subMod_gamm075.potentialSynapsesMotor_10))/(sum(block1_subMod_gamm075.potentialSynapsesMotor_10)+sum(block2_subMod_gamm075.potentialSynapsesInter_10))],...
                '-ko','LineWidth',2);
set(gca,'XTick',[1:4],'XTickLabels',[{'true'},2,5,10],'YLim',[1,10]);
ylabel('sum(off-Diag.)/sum(Diag)');
box off;
daspect([1,2,1])
%offsetAxes(gca);

subplot(4,4,4)
plot([1,2,3,4],[(sum(block1_subMod_gamm075.actualSynapsesMotor)+sum( block2_subMod_gamm075.actualSynapsesInter))/(sum(block1_subMod_gamm075.actualSynapsesInter)+sum(block2_subMod_gamm075.actualSynapsesMotor));...
                (sum(block1_subMod_gamm075.potentialSynapsesMotor_2)+sum(block2_subMod_gamm075.potentialSynapsesInter_2))/(sum(block2_subMod_gamm075.potentialSynapsesInter_2)+sum(block2_subMod_gamm075.potentialSynapsesMotor_2));...
                (sum(block1_subMod_gamm075.potentialSynapsesMotor_5)+sum(block2_subMod_gamm075.potentialSynapsesInter_5))/(sum(block2_subMod_gamm075.potentialSynapsesInter_5)+sum(block2_subMod_gamm075.potentialSynapsesMotor_5));...
                (sum(block1_subMod_gamm075.potentialSynapsesMotor_10)+sum(block2_subMod_gamm075.potentialSynapsesInter_10))/(sum(block2_subMod_gamm075.potentialSynapsesInter_10)+sum(block2_subMod_gamm075.potentialSynapsesMotor_10))],...
                '-ko','LineWidth',2);
set(gca,'XTick',[1:4],'XTickLabels',[{'true'},2,5,10]);
ylabel('sum(Diag.)/sum(off-Diag)');
box off;
daspect([10,1,1])

% figure;
% subplot(4,4,1)
% heatmap([mean(block1_subMod_gamm075.actualSynapsesMotor),mean(block1_subMod_gamm075.actualSynapsesInter);mean(block2_subMod_gamm075.actualSynapsesMotor),mean( block2_subMod_gamm075.actualSynapsesInter)],...
%     'ColorbarVisible','off','XData',{'ABDm','ABDi'},'YData',{'Int_ABDi','Int_ABDm'});
% title('True synapses');
% 
% subplot(4,4,2)
% heatmap([mean(block1_subMod_gamm075.potentialSynapsesMotor_5),mean(block1_subMod_gamm075.potentialSynapsesInter_5);mean(block2_subMod_gamm075.potentialSynapsesMotor_5),mean(block2_subMod_gamm075.potentialSynapsesInter_5)],...
%     'ColorbarVisible','off','XData',{'ABDm','ABDi'},'YData',{'Int_ABDi','Int_ABDm'});
% title('Potential synapses (5um)');


%% plot RS cell locations in block3, block4

Rov3 = [79961,81410,77449,83151,77456,77441,81611];
MiV2 = [78577,77694,77695,78566];
MiV1 = [80327,82221,82220,78940,82218,82217];
MiM1 = [77265];
RoM2 = [79244];
Mid3i = [76562];
RoL1 = [79395];

temp = distinguishable_colors(17,'w');
RSall = [Rov3,MiV2,MiV1,MiM1,RoM2,Mid3i,RoL1];
McellSoma = [66732, 25866, 16807];
McellSoma = TransformPoints(McellSoma,0);
transform_swc_AV(77099,temp(17,:),[],false,false);
transform_swc_AV(block3_gamma09_WithNull.cellIDs(ismember(block3_gamma09_WithNull.cellIDs,RSall)),temp(11,:),[],false,false);
transform_swc_AV(block4_gamma09_WithNull.cellIDs(ismember(block4_gamma09_WithNull.cellIDs,RSall)),temp(12,:),[],true,false);
scatter3(McellSoma(1),McellSoma(2),McellSoma(3),180,'MarkerEdgeColor',temp(17,:),'MarkerFaceColor',temp(17,:));


% transform_swc_AV(block3_gamma09_WithNull.cellIDs(ismember(block3_gamma09_WithNull.cellIDs,RSall)),hex2rgb('#91b4ba'),[],false,false);
% transform_swc_AV(block4_gamma09_WithNull.cellIDs(ismember(block4_gamma09_WithNull.cellIDs,RSall)),hex2rgb('#bb9a4e'),[],true,false);

%%

for i = 1:length(block2_gamma038_WithNull.cellIDs)
    A = SynapticPartners(block2_gamma038_WithNull.cellIDs(i),1,df);
    block2_gamma038_WithNull.recurrentFraction(i) = sum(ismember(A,block2_gamma038_WithNull.cellIDs))./length(A);
    block2_gamma038_WithNull.recurrentSynapses(i) = sum(ismember(A,block2_gamma038_WithNull.cellIDs));
    clear A;
end

for i = 1:length(block1_gamma038_WithNull.cellIDs)
    A = SynapticPartners(block1_gamma038_WithNull.cellIDs(i),1,df);
    block1_gamma038_WithNull.recurrentFraction(i) = sum(ismember(A,block1_gamma038_WithNull.cellIDs))./length(A);
    block1_gamma038_WithNull.recurrentSynapses(i)  = sum(ismember(A,block1_gamma038_WithNull.cellIDs));
    clear A;
end

% feedforward

figure;
subplot(4,4,1)
histogram(block2_gamma038_WithNull.recurrentFraction,'FaceColor',[0,0.5,1]);
hold on;
histogram(block1_gamma038_WithNull.recurrentFraction,'FaceColor',[1,0.5,0]);
legend({'Int','non-Int'});
axis square;
xlabel('% recurrent fraction');
box off;
offsetAxes(gca);
%%

[a,b,c] = OSI(block1_gamma038_WithNull.cellIDs,df);
[d,e,f] = OSI(block2_gamma038_WithNull.cellIDs,df);

g  = b+c;
h  = e+f;



block1_gamma038_WithNull.pdf =  fitdist(g(g>0),'Kernel');
block2_gamma038_WithNull.pdf = fitdist(h(h>0),'Kernel');

subplot(4,4,1)
% plot(0:1:100,pdf(block1_gamma038_WithNull.pdf,0:1:100),'-','color',[1,0.5,0],'LineWidth',2);
hold on;
% yyaxis right;
histogram(g(g>0),'BinWidth',5,'EdgeColor',[1,0.5,0],'DisplayStyle','stairs','LineWidth',2);
% yyaxis left
% plot(0:1:200,pdf(block2_gamma038_WithNull.pdf,0:1:200),'-','color',[0,0.5,1],'LineWidth',2);
% yyaxis right
histogram(h(h>0),'BinWidth',5,'EdgeColor',[0,0.5,1],'DisplayStyle','stairs','LineWidth',2);

axis square;
box off;
xlabel('ABD_M + ABD_I synapses');
ylabel('Count');
offsetAxes(gca);

[h1,p1] = kstest2(g(g>0),h(h>0));

%% Synapses size


% mean synapse size in block1_gamma038_WithNull (i-->j)

for i= 1:length(block1_gamma038_WithNull.cellIDs)
    if isExistReRoot(block1_gamma038_WithNull.cellIDs(i))
        block1_gamma038_WithNull.tree(i) = SwctoZbrian(block1_gamma038_WithNull.cellIDs(i));
    end
end

block1_gamma038_WithNull.b1tob1size = [];
block1_gamma038_WithNull.b1tob1postSynCoord = [];
block1_gamma038_WithNull.b1tob1postSynPathLength = [];
block1_gamma038_WithNull.b1tob1postSynPathLength_norm = [];

for i = 1:length(block1_gamma038_WithNull.cellIDs) % axon
    for j = 1:length(block1_gamma038_WithNull.cellIDs) % dendrite
        if ~isequal(block1_gamma038_WithNull.cellIDs(i),block1_gamma038_WithNull.cellIDs(j))
            temp = df.size(df.presyn_segid==block1_gamma038_WithNull.cellIDs(i) & df.postsyn_segid == block1_gamma038_WithNull.cellIDs(j));
            block1_gamma038_WithNull.b1tob1size = [block1_gamma038_WithNull.b1tob1size;temp];
            tempx = df.postsyn_x(df.presyn_segid==block1_gamma038_WithNull.cellIDs(i) & df.postsyn_segid == block1_gamma038_WithNull.cellIDs(j));
            tempy = df.postsyn_y(df.presyn_segid==block1_gamma038_WithNull.cellIDs(i) & df.postsyn_segid == block1_gamma038_WithNull.cellIDs(j));
            tempz = df.postsyn_z(df.presyn_segid==block1_gamma038_WithNull.cellIDs(i) & df.postsyn_segid == block1_gamma038_WithNull.cellIDs(j));
            block1_gamma038_WithNull.b1tob1postSynCoord  = [block1_gamma038_WithNull.b1tob1postSynCoord;[tempx,tempy,tempz]];
            if ~isempty(tempx)
                if isExistReRoot(block1_gamma038_WithNull.cellIDs(j))
                    block1_gamma038_WithNull.b1tob1postSynPathLength = [block1_gamma038_WithNull.b1tob1postSynPathLength;PathLengthToCoordinate(TransformPoints([tempx,tempy,tempz],0),block1_gamma038_WithNull.tree{j})];
                    normPathLength =  block1_gamma038_WithNull.b1tob1postSynPathLength ./max(Pvec_tree(block1_gamma038_WithNull.tree{j}));
                    block1_gamma038_WithNull.b1tob1postSynPathLength_norm = [block1_gamma038_WithNull.b1tob1postSynPathLength_norm;normPathLength];
                end
            end
            clear temp;
            clear tempx;
            clear tempy;
            clear tempz;
        end
    end
end


for i= 1:length(block2_gamma038_WithNull.cellIDs)
    if isExistReRoot(block2_gamma038_WithNull.cellIDs(i))
        block2_gamma038_WithNull.tree(i) = SwctoZbrian(block2_gamma038_WithNull.cellIDs(i));
    end
end

block1_gamma038_WithNull.b1tob2size = [];
block1_gamma038_WithNull.b1tob2postSynCoord = [];
block1_gamma038_WithNull.b1tob2postSynPathLength = [];
block1_gamma038_WithNull.b1tob2postSynPathLength_norm =[];

for i = 1:length(block1_gamma038_WithNull.cellIDs) % axon
    for j = 1:length(block2_gamma038_WithNull.cellIDs) % dendrite
        if ~isequal(block1_gamma038_WithNull.cellIDs(i),block2_gamma038_WithNull.cellIDs(j))
            temp = df.size(df.presyn_segid==block1_gamma038_WithNull.cellIDs(i) & df.postsyn_segid == block2_gamma038_WithNull.cellIDs(j));
            block1_gamma038_WithNull.b1tob2size  = [block1_gamma038_WithNull.b1tob2size ;temp];
            clear temp;
            tempx = df.postsyn_x(df.presyn_segid==block1_gamma038_WithNull.cellIDs(i) & df.postsyn_segid == block2_gamma038_WithNull.cellIDs(j));
            tempy = df.postsyn_y(df.presyn_segid==block1_gamma038_WithNull.cellIDs(i) & df.postsyn_segid == block2_gamma038_WithNull.cellIDs(j));
            tempz = df.postsyn_z(df.presyn_segid==block1_gamma038_WithNull.cellIDs(i) & df.postsyn_segid == block2_gamma038_WithNull.cellIDs(j));
            block1_gamma038_WithNull.b1tob2postSynCoord  = [block1_gamma038_WithNull.b1tob2postSynCoord;[tempx,tempy,tempz]];
            if ~isempty(tempx)
                if isExistReRoot(block2_gamma038_WithNull.cellIDs(j))
                    block1_gamma038_WithNull.b1tob2postSynPathLength = [block1_gamma038_WithNull.b1tob2postSynPathLength; PathLengthToCoordinate(TransformPoints([tempx,tempy,tempz],0),block2_gamma038_WithNull.tree{j})];
                    normPathLength =  block1_gamma038_WithNull.b1tob2postSynPathLength./max(Pvec_tree(block2_gamma038_WithNull.tree{j}));
                    block1_gamma038_WithNull.b1tob2postSynPathLength_norm = [block1_gamma038_WithNull.b1tob2postSynPathLength_norm;normPathLength];
                end
            end
            clear temp;
            clear tempx;
            clear tempy;
        end
    end
end

% mean synapse size in block2_gamma038_WithNull (i-->j)
block2_gamma038_WithNull.b2tob2size = [];
block2_gamma038_WithNull.b2tob2postSynCoord = [];
block2_gamma038_WithNull.b2tob2postSynPathLength = [];
block2_gamma038_WithNull.b2tob2postSynPathLength_norm=[];

for i = 1:length(block2_gamma038_WithNull.cellIDs) % axon
    for j = 1:length(block2_gamma038_WithNull.cellIDs) % dendrite
        if ~isequal(block2_gamma038_WithNull.cellIDs(i),block2_gamma038_WithNull.cellIDs(j))
            temp = df.size(df.presyn_segid==block2_gamma038_WithNull.cellIDs(i) & df.postsyn_segid == block2_gamma038_WithNull.cellIDs(j));
            block2_gamma038_WithNull.b2tob2size = [block2_gamma038_WithNull.b2tob2size;temp];
            clear temp;
            tempx = df.postsyn_x(df.presyn_segid==block2_gamma038_WithNull.cellIDs(i) & df.postsyn_segid == block2_gamma038_WithNull.cellIDs(j));
            tempy = df.postsyn_y(df.presyn_segid==block2_gamma038_WithNull.cellIDs(i) & df.postsyn_segid == block2_gamma038_WithNull.cellIDs(j));
            tempz = df.postsyn_z(df.presyn_segid==block2_gamma038_WithNull.cellIDs(i) & df.postsyn_segid == block2_gamma038_WithNull.cellIDs(j));
            block2_gamma038_WithNull.b2tob2postSynCoord  = [block2_gamma038_WithNull.b2tob2postSynCoord;[tempx,tempy,tempz]];
            if ~isempty(tempx)
                if isExistReRoot(block2_gamma038_WithNull.cellIDs(j))
                    block2_gamma038_WithNull.b2tob2postSynPathLength = [block2_gamma038_WithNull.b2tob2postSynPathLength;PathLengthToCoordinate(TransformPoints([tempx,tempy,tempz],0),block2_gamma038_WithNull.tree{j})];
                    normPathLength =   block2_gamma038_WithNull.b2tob2postSynPathLength./max(Pvec_tree(block2_gamma038_WithNull.tree{j}));
                    block2_gamma038_WithNull.b2tob2postSynPathLength_norm = [block2_gamma038_WithNull.b2tob2postSynPathLength_norm;normPathLength];
                end
            end
            clear temp;
            clear tempx;
            clear tempy;
        end
    end
end


block2_gamma038_WithNull.b2tob1size = [];
block2_gamma038_WithNull.b2tob1postSynCoord = [];
block2_gamma038_WithNull.b2tob1postSynPathLength =[];
block2_gamma038_WithNull.b2tob1postSynPathLength_norm = [];


for i = 1:length(block2_gamma038_WithNull.cellIDs) % axon
    for j = 1:length(block1_gamma038_WithNull.cellIDs) % dendrite
        if ~isequal(block2_gamma038_WithNull.cellIDs(i),block1_gamma038_WithNull.cellIDs(j))
            temp = df.size(df.presyn_segid==block2_gamma038_WithNull.cellIDs(i) & df.postsyn_segid == block1_gamma038_WithNull.cellIDs(j));
            block2_gamma038_WithNull.b2tob1size = [block2_gamma038_WithNull.b2tob1size;temp];
            clear temp;
            tempx = df.postsyn_x(df.presyn_segid==block2_gamma038_WithNull.cellIDs(i) & df.postsyn_segid == block1_gamma038_WithNull.cellIDs(j));
            tempy = df.postsyn_y(df.presyn_segid==block2_gamma038_WithNull.cellIDs(i) & df.postsyn_segid == block1_gamma038_WithNull.cellIDs(j));
            tempz = df.postsyn_z(df.presyn_segid==block2_gamma038_WithNull.cellIDs(i) & df.postsyn_segid == block1_gamma038_WithNull.cellIDs(j));
            block2_gamma038_WithNull.b2tob1postSynCoord  = [block2_gamma038_WithNull.b2tob1postSynCoord;[tempx,tempy,tempz]];
            if ~isempty(tempx)
                if isExistReRoot(block1_gamma038_WithNull.cellIDs(j))
                    block2_gamma038_WithNull.b2tob1postSynPathLength = [block2_gamma038_WithNull.b2tob1postSynPathLength;PathLengthToCoordinate(TransformPoints([tempx,tempy,tempz],0),block1_gamma038_WithNull.tree{j})];
                    normPathLength =   block2_gamma038_WithNull.b2tob1postSynPathLength./max(Pvec_tree(block1_gamma038_WithNull.tree{j}));
                    block2_gamma038_WithNull.b2tob1postSynPathLength_norm = [block2_gamma038_WithNull.b2tob1postSynPathLength_norm;normPathLength];
                end
            end
            clear tempz;
            clear tempx;
            clear tempy;
        end
    end
end

%%
figure;
subplot(4,4,1)
histogram(log10(block1_gamma038_WithNull.b1tob1size),'BinWidth',0.1,'EdgeColor',[1,0.5,0],'DisplayStyle','stairs','LineWidth',2);
hold on;
histogram(log10(block2_gamma038_WithNull.b2tob2size),'BinWidth',0.1,'EdgeColor',[0,0.5,1],'DisplayStyle','stairs','LineWidth',2);
axis square;
box off;
set(gca,'XTick',[2,3,4],'XTickLabel',[10^2,10^3,10^4]);
xlabel('log size (vx)');
ylabel('count');
legend({'mod1->mod1','mod2-mod2'},'Location','bestoutside');
% line([mean(log10(block1_gamma038_WithNull.b1tob1size)),mean(log10(block1_gamma038_WithNull.b1tob1size))],...
%     [0,500],'color',[1,0.5,0],'lineWidth',2);
% line([mean(log10(block2_gamma038_WithNull.b2tob2size)),mean(log10(block2_gamma038_WithNull.b2tob2size))],...
%     [0,500],'color',[0,0.5,1],'lineWidth',2);

subplot(4,4,2)
histogram(log10(block1_gamma038_WithNull.b1tob2size),'BinWidth',0.1,'EdgeColor',[1,0.5,0],'DisplayStyle','stairs','LineWidth',2,'EdgeAlpha',0.4);
hold on;
histogram(log10(block2_gamma038_WithNull.b2tob1size),'BinWidth',0.1,'EdgeColor',[0,0.5,1],'DisplayStyle','stairs','LineWidth',2,'EdgeAlpha',0.4);
axis square;
box off;
set(gca,'XTick',[2,3,4],'XTickLabel',[10^2,10^3,10^4]);
legend({'mod1->mod2','mod2->mod1'},'Location','bestoutside');
%%

% mean synapse size in block2_gamma038_WithNull to motor (i-->j)
ABDr_CellIDs = [82140,82145,82143,77648,82146,77302,77301,82194,77705,77710,77305,77672,77300,77709,77661];
ABDc_CellIDs = [82213,77654,77646,82212,77295,81172,77657,82197,77682,77154, 77652,77658,77628,77292,77688,82195,77296];
ABDIr_CellIDs = [78553,77150,77886,78547,77631,77158,77665,78552,77668,77618,77634,78556];
ABDIc_CellIDs = [78574,79051,77148,77625,77144,77692,77641,79066,77640,77643];

ABDmotor = [ABDr_CellIDs,ABDc_CellIDs];
ABDinter  = [ABDIr_CellIDs,ABDIc_CellIDs];

block2_gamma038_WithNull.b2toABDm = [];
block2_gamma038_WithNull.b2toABDmPathLength = [];
block2_gamma038_WithNull.b2toABDmPostSynCoord = [];

for i = 1:length(block2_gamma038_WithNull.cellIDs) % axon
    for j = 1:length(ABDmotor) % dendrite
        if ~isequal(block2_gamma038_WithNull.cellIDs(i),ABDmotor(j))
            temp = df.size(df.presyn_segid==block2_gamma038_WithNull.cellIDs(i) & df.postsyn_segid == ABDmotor(j));
            block2_gamma038_WithNull.b2toABDm = [block2_gamma038_WithNull.b2toABDm;temp];
            clear temp;
            tempx = df.postsyn_x(df.presyn_segid==block2_gamma038_WithNull.cellIDs(i) & df.postsyn_segid == ABDmotor(j));
            tempy = df.postsyn_y(df.presyn_segid==block2_gamma038_WithNull.cellIDs(i) & df.postsyn_segid == ABDmotor(j));
            tempz = df.postsyn_z(df.presyn_segid==block2_gamma038_WithNull.cellIDs(i) & df.postsyn_segid == ABDmotor(j));
            block2_gamma038_WithNull.b2toABDmPostSynCoord  = [block2_gamma038_WithNull.b2toABDmPostSynCoord;[tempx,tempy,tempz]];
            if ~isempty(tempx)
                if isExistReRoot(block1_gamma038_WithNull.cellIDs(j))
                    block2_gamma038_WithNull.b2toABDmPathLength = [block2_gamma038_WithNull.b2toABDmPathLength;PathLengthToCoordinate(TransformPoints([tempx,tempy,tempz],0),block1_gamma038_WithNull.tree{j})];
%                     normPathLength =   block2_gamma038_WithNull.b2tob1postSynPathLength./max(Pvec_tree(block1_gamma038_WithNull.tree{j}));
%                     block2_gamma038_WithNull.b2tob1postSynPathLength_norm = [block2_gamma038_WithNull.b2tob1postSynPathLength_norm;normPathLength];
                end
            end
            clear tempz;
            clear tempx;
            clear tempy;
            
        end
    end
end

block2_gamma038_WithNull.b2toABDi = [];
block2_gamma038_WithNull.b2toABDiPathLength = [];
block2_gamma038_WithNull.b2toABDiPostSynCoord = [];

for i = 1:length(block2_gamma038_WithNull.cellIDs) % axon
    for j = 1:length(ABDinter) % dendrite
        if ~isequal(block2_gamma038_WithNull.cellIDs(i),ABDinter(j))
            temp = df.size(df.presyn_segid==block2_gamma038_WithNull.cellIDs(i) & df.postsyn_segid == ABDinter(j));
            block2_gamma038_WithNull.b2toABDi = [block2_gamma038_WithNull.b2toABDi;temp];
            clear temp;
            tempx = df.postsyn_x(df.presyn_segid==block2_gamma038_WithNull.cellIDs(i) & df.postsyn_segid == ABDinter(j));
            tempy = df.postsyn_y(df.presyn_segid==block2_gamma038_WithNull.cellIDs(i) & df.postsyn_segid == ABDinter(j));
            tempz = df.postsyn_z(df.presyn_segid==block2_gamma038_WithNull.cellIDs(i) & df.postsyn_segid == ABDinter(j));
            block2_gamma038_WithNull.b2toABDiPostSynCoord  = [block2_gamma038_WithNull.b2toABDiPostSynCoord;[tempx,tempy,tempz]];
            if ~isempty(tempx)
                if isExistReRoot(block1_gamma038_WithNull.cellIDs(j))
                    block2_gamma038_WithNull.b2toABDiPathLength = [block2_gamma038_WithNull.b2toABDiPathLength;PathLengthToCoordinate(TransformPoints([tempx,tempy,tempz],0),block1_gamma038_WithNull.tree{j})];
%                     normPathLength =   block2_gamma038_WithNull.b2tob1postSynPathLength./max(Pvec_tree(block1_gamma038_WithNull.tree{j}));
%                     block2_gamma038_WithNull.b2tob1postSynPathLength_norm = [block2_gamma038_WithNull.b2tob1postSynPathLength_norm;normPathLength];
                end
            end
            clear tempz;
            clear tempx;
            clear tempy;
        end
    end
end

%%
figure;

subplot(4,4,1)
histogram(log10(block2_gamma038_WithNull.b2toABDm),'BinWidth',0.1,'EdgeColor',block2_subMod_gamm075.color,'DisplayStyle','stairs','LineWidth',2);
hold on;
histogram(log10(block2_gamma038_WithNull.b2toABDi),'BinWidth',0.1,'EdgeColor',block1_subMod_gamm075.color,'DisplayStyle','stairs','LineWidth',2);
box off;
axis square;
set(gca,'XTick',[2,3,4],'XTickLabel',[10^2,10^3,10^4]);
legend({'mod2b->ABD_m','mod2a-> ABD_i'},'Location','bestoutside');
xlabel('log size (vx)');
ylabel('count');


subplot(4,4,2)
histogram(log10(block2_gamma038_WithNull.b2tob2size),'BinWidth',0.1,'EdgeColor',[0,0.5,1],'DisplayStyle','stairs','LineWidth',2);

hold on;
histogram(log10(block2_gamma038_WithNull.b2toABDm),'BinWidth',0.1,'EdgeColor','g','DisplayStyle','stairs','LineWidth',2);
histogram(log10(block2_gamma038_WithNull.b2toABDi),'BinWidth',0.1,'EdgeColor','m','DisplayStyle','stairs','LineWidth',2);
box off;
axis square;
set(gca,'XTick',[2,3,4],'XTickLabel',[10^2,10^3,10^4]);
legend({'b2->b2','b2->ABD_m','b2-> ABD_i'},'Location','bestoutside');


%%

subplot(4,4,1)
histogram(block1_gamma038_WithNull.b1tob1postSynPathLength,'BinWidth',10,'EdgeColor','none','FaceColor',[1,0.5,0],'Normalization','probability');
hold on
histogram(block2_gamma038_WithNull.b2tob1postSynPathLength,'BinWidth',10,'EdgeColor',[0,0.5,1],'FaceColor','none','LineWidth',2,'DisplayStyle','stairs','Normalization','probability');

axis square;
box off;
xlabel('Pathlength (\mum)');
ylabel('count');
legend({'mod1->mod1','mod2->mod1'},'Location','bestoutside');


subplot(4,4,2)
histogram(block2_gamma038_WithNull.b2tob2postSynPathLength,'BinWidth',10,'EdgeColor','none','FaceColor',[0,0.5,1],'Normalization','probability');
hold on;
histogram(block1_gamma038_WithNull.b1tob2postSynPathLength,'BinWidth',10,'EdgeColor',[1,0.5,0],'FaceColor','none','LineWidth',2,'DisplayStyle','stairs','Normalization','probability');

axis square;
box off;
xlabel('Pathlength (\mum)');
ylabel('count');
legend({'mod2->mod2','mod1->mod2'},'Location','bestoutside');

subplot(4,4,5)
histogram(block2_gamma038_WithNull.b2toABDmPathLength,'BinWidth',12,'EdgeColor','none','EdgeColor',block2_subMod_gamm075.color,...
    'Normalization','probability','DisplayStyle','stairs','LineWidth',2);
hold on
histogram(block2_gamma038_WithNull.b2toABDiPathLength,'BinWidth',12,'EdgeColor','none','EdgeColor',block1_subMod_gamm075.color,...
    'Normalization','probability','DisplayStyle','stairs','LineWidth',2);
axis square;
box off;
xlabel('Pathlength (\mum)');
ylabel('count');

%% dual synapses coorelations
block1_gamma038_WithNull.pairs_all = [];
block1_gamma038_WithNull.pairs_single = [];
block1_gamma038_WithNull.pairs_dual = [];
block1_gamma038_WithNull.pairs_triple = [];
block1_gamma038_WithNull.pairs_quart = [];
block1_gamma038_WithNull.pairs_pent = [];

for i = 1:length(block1_gamma038_WithNull.cellIDs)
    for j = 1:length(block1_gamma038_WithNull.cellIDs)
        if ~isequal(block1_gamma038_WithNull.cellIDs(i),block1_gamma038_WithNull.cellIDs(j))
            sz = df.size(df.presyn_segid==block1_gamma038_WithNull.cellIDs(i) & df.postsyn_segid == block1_gamma038_WithNull.cellIDs(j));
            block1_gamma038_WithNull.pairs_all = [block1_gamma038_WithNull.pairs_all;repmat(block1_gamma038_WithNull.cellIDs(i),length(sz),1),repmat(block1_gamma038_WithNull.cellIDs(j),length(sz),1),sz];
            if size(sz,1) == 1
                block1_gamma038_WithNull.pairs_single =  [block1_gamma038_WithNull.pairs_single;repmat(block1_gamma038_WithNull.cellIDs(i),length(sz),1),repmat(block1_gamma038_WithNull.cellIDs(j),length(sz),1),sz];
            elseif size(sz,1) == 2
                block1_gamma038_WithNull.pairs_dual =  [block1_gamma038_WithNull.pairs_dual;repmat(block1_gamma038_WithNull.cellIDs(i),length(sz),1),repmat(block1_gamma038_WithNull.cellIDs(j),length(sz),1),sz];
            elseif size(sz,1) == 3
                block1_gamma038_WithNull.pairs_triple =  [block1_gamma038_WithNull.pairs_triple;repmat(block1_gamma038_WithNull.cellIDs(i),length(sz),1),repmat(block1_gamma038_WithNull.cellIDs(j),length(sz),1),sz];
            elseif size(sz,1) == 4
                block1_gamma038_WithNull.pairs_quart =  [block1_gamma038_WithNull.pairs_quart;repmat(block1_gamma038_WithNull.cellIDs(i),length(sz),1),repmat(block1_gamma038_WithNull.cellIDs(j),length(sz),1),sz];
            else
                block1_gamma038_WithNull.pairs_pent =  [block1_gamma038_WithNull.pairs_pent;repmat(block1_gamma038_WithNull.cellIDs(i),length(sz),1),repmat(block1_gamma038_WithNull.cellIDs(j),length(sz),1),sz];
            end
            clear sz
        end
    end
end

%scatter(block1_gamma038_WithNull.pairs_dual(1:2:end,3),block1_gamma038_WithNull.pairs_dual(2:2:end,3))
block1_gamma038_WithNull.dualpairs = [];
for i = 1:length(block1_gamma038_WithNull.pairs_dual)
    block1_gamma038_WithNull.dualpairs = [block1_gamma038_WithNull.dualpairs;nchoosek(block1_gamma038_WithNull.pairs_dual(i:i+1,3),2)];
end
block1_gamma038_WithNull.dualtriple = [];
for i = 1:length(block1_gamma038_WithNull.pairs_triple)
    block1_gamma038_WithNull.dualtriple = [block1_gamma038_WithNull.dualtriple;nchoosek(block1_gamma038_WithNull.pairs_triple(i:i+2,3),2)];
end
block1_gamma038_WithNull.dualquart = [];
for i = 1:length(block1_gamma038_WithNull.pairs_quart)
    block1_gamma038_WithNull.dualquart = [block1_gamma038_WithNull.dualquart;nchoosek(block1_gamma038_WithNull.pairs_quart(i:i+3,3),2)];
end
block1_gamma038_WithNull.dualpent = [];
for i = 1:length(block1_gamma038_WithNull.pairs_pent)
    block1_gamma038_WithNull.dualpent = [block1_gamma038_WithNull.dualpent;nchoosek(block1_gamma038_WithNull.pairs_pent(i:i+4,3),2)];
end
% 
% block1_gamma038_WithNull.diffPairs = abs(block1_gamma038_WithNull.pairs_all(:,1)-block1_gamma038_WithNull.pairs_all(:,2))./block1_gamma038_WithNull.pairs_all(:,1); % all pairs
% block1_gamma038_WithNull.DualdiffPairs = abs(block1_gamma038_WithNull.pairs_dual(:,1)-block1_gamma038_WithNull.pairs_dual(:,2)); % only dual pairs
% block1_gamma038_WithNull.locs = diff(block1_gamma038_WithNull.DualdiffPairs);
% block1_gamma038_WithNull.index = find(block1_gamma038_WithNull.locs);

% ind = 1;
% block1_gamma038_WithNull.dualPairs = [];
% for i = 1:length(block1_gamma038_WithNull.index)
%     nums = block1_gamma038_WithNull.index(i)-ind;
%     block1_gamma038_WithNull.dualPairs = [block1_gamma038_WithNull.DualdiffPairs;nchoosek(block1_gamma038_WithNull.pairs_dual(ind:block1_gamma038_WithNull.index(i),3),2)];
%     ind = block1_gamma038_WithNull.index(i)+1;
% end


% block1_gamma038_WithNull-->b2
block1_gamma038_WithNull.b1tob2pairs_all = [];
block1_gamma038_WithNull.b1tob2pairs_dual = [];

for i = 1:length(block1_gamma038_WithNull.cellIDs)
    for j = 1:length(block2_gamma038_WithNull.cellIDs)
        if ~isequal(block1_gamma038_WithNull.cellIDs(i),block2_gamma038_WithNull.cellIDs(j))
            sz = df.size(df.presyn_segid==block1_gamma038_WithNull.cellIDs(i) & df.postsyn_segid == block2_gamma038_WithNull.cellIDs(j));
            block1_gamma038_WithNull.b1tob2pairs_all = [block1_gamma038_WithNull.b1tob2pairs_all;repmat(block1_gamma038_WithNull.cellIDs(i),length(sz),1),repmat(block2_gamma038_WithNull.cellIDs(j),length(sz),1),sz];
            if size(sz,1)>1
                block1_gamma038_WithNull.b1tob2pairs_dual =  [block1_gamma038_WithNull.b1tob2pairs_dual;repmat(block1_gamma038_WithNull.cellIDs(i),length(sz),1),repmat(block2_gamma038_WithNull.cellIDs(j),length(sz),1),sz];
            end
            clear sz
        end
    end
end

block1_gamma038_WithNull.b1tob2diffPairs = abs(block1_gamma038_WithNull.b1tob2pairs_all(:,1)-block1_gamma038_WithNull.b1tob2pairs_all(:,2))./block1_gamma038_WithNull.b1tob2pairs_all(:,1);



% % gmm for dual pairs
% 
% % monovariate
% block1_gamma038_WithNull.gmmDualpairs = fitgmdist(log10(block1_gamma038_WithNull.pairs_dual(:,3)),2,'Options',statset('MaxIter',1000));
% samples = 2:0.01:5;
% evalat = samples(randi(length(samples),1000,1))';
% block1_gamma038_WithNull.gmmDualpairsPDF = pdf(block1_gamma038_WithNull.gmmDualpairs,[2:0.01:5]');
% 
% %bivariate
% 
% block1_gamma038_WithNull.gmmDualpairsBI = fitgmdist([log10(block1_gamma038_WithNull.dualPairs(:,1)),log10(block1_gamma038_WithNull.dualPairs(:,2))],2,'Options',statset('MaxIter',1000));
% block1_gamma038_WithNull.gmmDualpairsBIPDF = @(x,y)reshape(pdf(block1_gamma038_WithNull.gmmDualpairsBI,[x(:),y(:)]),size(x));
% 
% 
% subplot(2,2,1);
% histogram(log10(block1_gamma038_WithNull.pairs_dual(:,3)),'BinWidth',0.1,'EdgeColor',[1,0.5,0],'DisplayStyle','stairs','LineWidth',2);
% xlabel('size(vx)');
% ylabel('count');
% set(gca,'XTick',[2,3,4,5],'XTickLabels',[10^2,10^3,10^4,10^5]);
% hold on
% yyaxis right;
% plot(2:0.01:5,block1_gamma038_WithNull.gmmDualpairsPDF,'color',[1,0.5,0],'LineWidth',2);
% axis square;
% 
% subplot(2,2,2)
% scatter(log10(block1_gamma038_WithNull.dualPairs(:,1)),log10(block1_gamma038_WithNull.dualPairs(:,2)),'.','markerEdgecolor',[1,0.5,0],'MarkerEdgeAlpha',0.2);
% hold on;
% fcontour(block1_gamma038_WithNull.gmmDualpairsBIPDF,[2 5],'LineWidth',2);
% set(gca,'XTick',[2,3,4,5],'XTickLabels',[10^2,10^3,10^4,10^5],'YTick',[2,3,4,5],'YTickLabels',[10^2,10^3,10^4,10^5]);
% 
% axis square;

block2_gamma038_WithNull.pairs_all = [];
block2_gamma038_WithNull.pairs_single = [];
block2_gamma038_WithNull.pairs_dual = [];
block2_gamma038_WithNull.pairs_triple = [];
block2_gamma038_WithNull.pairs_quart = [];
block2_gamma038_WithNull.pairs_pent = [];

for i = 1:length(block2_gamma038_WithNull.cellIDs)
    for j = 1:length(block2_gamma038_WithNull.cellIDs)
        if ~isequal(block2_gamma038_WithNull.cellIDs(i),block2_gamma038_WithNull.cellIDs(j))
            sz = df.size(df.presyn_segid==block2_gamma038_WithNull.cellIDs(i) & df.postsyn_segid == block2_gamma038_WithNull.cellIDs(j));
            block2_gamma038_WithNull.pairs_all = [block2_gamma038_WithNull.pairs_all;repmat(block2_gamma038_WithNull.cellIDs(i),length(sz),1),repmat(block2_gamma038_WithNull.cellIDs(j),length(sz),1),sz];
            if size(sz,1) == 1
                block2_gamma038_WithNull.pairs_single =  [block2_gamma038_WithNull.pairs_single;repmat(block2_gamma038_WithNull.cellIDs(i),length(sz),1),repmat(block2_gamma038_WithNull.cellIDs(j),length(sz),1),sz];
            elseif size(sz,1) == 2
                block2_gamma038_WithNull.pairs_dual =  [block2_gamma038_WithNull.pairs_dual;repmat(block2_gamma038_WithNull.cellIDs(i),length(sz),1),repmat(block2_gamma038_WithNull.cellIDs(j),length(sz),1),sz];
            elseif size(sz,1) == 3
                block2_gamma038_WithNull.pairs_triple =  [block2_gamma038_WithNull.pairs_triple;repmat(block2_gamma038_WithNull.cellIDs(i),length(sz),1),repmat(block2_gamma038_WithNull.cellIDs(j),length(sz),1),sz];
            elseif size(sz,1) == 4
                block2_gamma038_WithNull.pairs_quart =  [block2_gamma038_WithNull.pairs_quart;repmat(block2_gamma038_WithNull.cellIDs(i),length(sz),1),repmat(block2_gamma038_WithNull.cellIDs(j),length(sz),1),sz];
            else
                block2_gamma038_WithNull.pairs_pent =  [block2_gamma038_WithNull.pairs_pent;repmat(block2_gamma038_WithNull.cellIDs(i),length(sz),1),repmat(block2_gamma038_WithNull.cellIDs(j),length(sz),1),sz];
            end
            clear sz
        end
    end
end

block2_gamma038_WithNull.dualpairs = [];
for i = 1:length(block1_gamma038_WithNull.pairs_dual)
    block2_gamma038_WithNull.dualpairs = [block2_gamma038_WithNull.dualpairs;nchoosek(block2_gamma038_WithNull.pairs_dual(i:i+1,3),2)];
end
block2_gamma038_WithNull.dualtriple = [];
for i = 1:length(block2_gamma038_WithNull.pairs_triple)
    block2_gamma038_WithNull.dualtriple = [block2_gamma038_WithNull.dualtriple;nchoosek(block2_gamma038_WithNull.pairs_triple(i:i+2,3),2)];
end
block2_gamma038_WithNull.dualquart = [];
for i = 1:length(block2_gamma038_WithNull.pairs_quart)
    block2_gamma038_WithNull.dualquart = [block2_gamma038_WithNull.dualquart;nchoosek(block2_gamma038_WithNull.pairs_quart(i:i+3,3),2)];
end
block2_gamma038_WithNull.dualpent = [];
for i = 1:length(block2_gamma038_WithNull.pairs_pent)
    block2_gamma038_WithNull.dualpent = [block2_gamma038_WithNull.dualpent;nchoosek(block2_gamma038_WithNull.pairs_pent(i:i+4,3),2)];
end


% block2_gamma038_WithNull.diffPairs = abs(block2_gamma038_WithNull.pairs_all(:,1)-block2_gamma038_WithNull.pairs_all(:,2));
% block2_gamma038_WithNull.DualdiffPairs = abs(block2_gamma038_WithNull.pairs_dual(:,1)-block2_gamma038_WithNull.pairs_dual(:,2));
% 
% block2_gamma038_WithNull.locs = diff(block2_gamma038_WithNull.diffPairs);
% block2_gamma038_WithNull.index = find(block2_gamma038_WithNull.locs);
% 
% ind = 1;
% block2_gamma038_WithNull.dualPairs = [];
% for i = 1:length(block2_gamma038_WithNull.index)
%     nums = block2_gamma038_WithNull.index(i)-ind;
%     block2_gamma038_WithNull.dualPairs = [block2_gamma038_WithNull.dualPairs;nchoosek(block2_gamma038_WithNull.pairs_dual(ind:block2_gamma038_WithNull.index(i),3),2)];
%     ind = block2_gamma038_WithNull.index(i)+1;
% end
% 
%b2-->b1
block2_gamma038_WithNull.b2tob1pairs_all = [];
block2_gamma038_WithNull.b2tob1pairs_dual = [];

for i = 1:length(block2_gamma038_WithNull.cellIDs)
    for j = 1:length(block1_gamma038_WithNull.cellIDs)
        if ~isequal(block2_gamma038_WithNull.cellIDs(i),block1_gamma038_WithNull.cellIDs(j))
            sz = df.size(df.presyn_segid==block2_gamma038_WithNull.cellIDs(i) & df.postsyn_segid == block1_gamma038_WithNull.cellIDs(j));
            block2_gamma038_WithNull.b2tob1pairs_all = [block2_gamma038_WithNull.b2tob1pairs_all;repmat(block2_gamma038_WithNull.cellIDs(i),length(sz),1),repmat(block1_gamma038_WithNull.cellIDs(j),length(sz),1),sz];
            if size(sz,1)>1
                block2_gamma038_WithNull.b2tob1pairs_dual =  [block2_gamma038_WithNull.b2tob1pairs_dual;repmat(block2_gamma038_WithNull.cellIDs(i),length(sz),1),repmat(block1_gamma038_WithNull.cellIDs(j),length(sz),1),sz];
            end
            clear sz
        end
    end
end

block2_gamma038_WithNull.b2tob1diffPairs = abs(block2_gamma038_WithNull.b2tob1pairs_all(:,1)-block2_gamma038_WithNull.b2tob1pairs_all(:,2));

%block2_gamma038_WithNull.DualdiffPairs = abs(block2_gamma038_WithNull.pairs_dual(:,1)-block2_gamma038_WithNull.pairs_dual(:,2));



figure;
subplot(4,4,1)
scatter(log10(block1_gamma038_WithNull.dualPairs(:,1)),log10(block1_gamma038_WithNull.dualPairs(:,2)),'.','MarkerEdgeColor',[1,0.5,0],'MarkerEdgeAlpha',0.2);
hold on
f = showfit(ezfit(log10(block1_gamma038_WithNull.dualPairs(:,1)),log10(block1_gamma038_WithNull.dualPairs(:,2)),'affine'),'fitcolor',[1,0.5,0],'dispeqboxmode','off');
text(2.5,4,sprintf('r = %0.4f',f.r));
axis square;
set(gca,'XTick',[2,3,4,5],'XTickLabels',[10^2,10^3,10^4,10^5],'YTick',[2,3,4,5],'YTickLabels',[10^2,10^3,10^4,10^5]);

subplot(4,4,2)
scatter(log10(block2_gamma038_WithNull.dualPairs(:,1)),log10(block2_gamma038_WithNull.dualPairs(:,2)),'.','MarkerEdgeColor',[0,0.5,1],'MarkerEdgeAlpha',0.2);
hold on
f = showfit(ezfit(log10(block2_gamma038_WithNull.dualPairs(:,1)),log10(block2_gamma038_WithNull.dualPairs(:,2)),'affine'),'fitcolor',[0,0.5,1],'dispeqboxmode','off');
text(2.5,4,sprintf('r = %0.4f',f.r));
axis square;
set(gca,'XTick',[2,3,4,5],'XTickLabels',[10^2,10^3,10^4,10^5],'YTick',[2,3,4,5],'YTickLabels',[10^2,10^3,10^4,10^5]);

% h = scatterhist([log10(block1_gamma038_WithNull.dualPairs(:,1));log10(block2_gamma038_WithNull.dualPairs(:,1))],...
%     [log10(block1_gamma038_WithNull.dualPairs(:,2));log10(block2_gamma038_WithNull.dualPairs(:,2))],...
%     'Group',[ones(size(block1_gamma038_WithNull.dualPairs,1),1);2*ones(size(block2_gamma038_WithNull.dualPairs,1),1)],'Kernel','on','color',[[1,0.5,0];[0,0.5,1]],...
%     'marker','.');
% set(h(1), 'XTick',[2,3,4,5],'XTickLabel',[10^2,10^3,10^4],'YTick',[2,3,4,5],'YTickLabel',[10^2,10^3,10^4]);
figure;
subplot(4,4,1)
histogram(histcounts(block1_gamma038_WithNull.b1tob2diffPairs,unique( block1_gamma038_WithNull.b1tob2diffPairs)),'Normalization','cdf','DisplayStyle','stairs','EdgeColor',[1,0.5,0],'LineWidth',2);
hold on; histogram(histcounts( block1_gamma038_WithNull.diffPairs,unique( block1_gamma038_WithNull.diffPairs)),'Normalization','cdf','FaceColor',[1,0.5,0],'LineWidth',2,'EdgeColor','none');
axis square;
xlabel('number of synaptic pairs');
ylabel('probability');
legend({'b1->b2','b1->b1'});

subplot(4,4,2)
histogram(histcounts(block2_gamma038_WithNull.b2tob1diffPairs,unique( block2_gamma038_WithNull.b2tob1diffPairs)),'Normalization','cdf','DisplayStyle','stairs','EdgeColor',[0,0.5,1],'LineWidth',2);
hold on; histogram(histcounts( block2_gamma038_WithNull.diffPairs,unique( block2_gamma038_WithNull.diffPairs)),'Normalization','cdf','FaceColor',[0,0.5,1],'LineWidth',2,'EdgeColor','none');
axis square;
legend({'b2->b1','b2->b2'});
%%

figure

% subplot(2,5,1)
% scatter(block1_gamma038_WithNull.pairs_single(:,1),block1_gamma038_WithNull.pairs_single(:,2),'MarkerFaceColor',[1,0.5,0],'MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
% hold on;
% f = showfit(ezfit(block1_gamma038_WithNull.pairs_single(:,1),block1_gamma038_WithNull.pairs_single(:,2),'affine'),'fitcolor','k','dispeqboxmode','off');
% %set(gca, 'XScale','log','YScale','log');
% axis square;
% text(10e2,10e3,sprintf('r = %0.2f',f.r));
% 

subplot(2,5,2)
scatter(block1_gamma038_WithNull.dualpairs(:,1),block1_gamma038_WithNull.dualpairs(:,2),'MarkerFaceColor',[1,0.5,0],'MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
hold on;
f = showfit(ezfit(block1_gamma038_WithNull.dualpairs(:,1),block1_gamma038_WithNull.dualpairs(:,2),'affine'),'fitcolor','k','dispeqboxmode','off');
%set(gca, 'XScale','log','YScale','log');
axis square;
text(10e2,10e3,sprintf('r = %0.2f',f.r));

subplot(2,5,3)
scatter(block1_gamma038_WithNull.dualtriple(:,1),block1_gamma038_WithNull.dualtriple(:,2),'MarkerFaceColor',[1,0.5,0],'MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
hold on;
f = showfit(ezfit(block1_gamma038_WithNull.dualtriple(:,1),block1_gamma038_WithNull.dualtriple(:,2),'affine'),'fitcolor','k','dispeqboxmode','off');
%set(gca, 'XScale','log','YScale','log');
axis square;
text(10e2,10e3,sprintf('r = %0.2f',f.r));

subplot(2,5,4)
scatter(block1_gamma038_WithNull.dualquart(:,1),block1_gamma038_WithNull.dualquart(:,2),'MarkerFaceColor',[1,0.5,0],'MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
hold on;
f = showfit(ezfit(block1_gamma038_WithNull.dualquart(:,1),block1_gamma038_WithNull.dualquart(:,2),'affine'),'fitcolor','k','dispeqboxmode','off');
%set(gca, 'XScale','log','YScale','log');
axis square;
text(10e2,10e3,sprintf('r = %0.2f',f.r));

subplot(2,5,5)
scatter(block1_gamma038_WithNull.dualpent(:,1),block1_gamma038_WithNull.dualpent(:,2),'MarkerFaceColor',[1,0.5,0],'MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
hold on;
f = showfit(ezfit(block1_gamma038_WithNull.dualpent(:,1),block1_gamma038_WithNull.dualpent(:,2),'affine'),'fitcolor','k','dispeqboxmode','off');
%set(gca, 'XScale','log','YScale','log');
axis square;
text(10e2,10e3,sprintf('r = %0.2f',f.r));



% subplot(2,5,6)
% scatter(block2_gamma038_WithNull.dualpairs(:,1),block2_gamma038_WithNull.dualpairs(:,2),'MarkerFaceColor',[0,0.5,1],'MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
% hold on;
% f = showfit(ezfit(block2_gamma038_WithNull.dualpairs(:,1),block2_gamma038_WithNull.dualpairs(:,2),'affine'),'fitcolor','k','dispeqboxmode','off');
% %set(gca, 'XScale','log','YScale','log');
% axis square;
% text(10e2,10e3,sprintf('r = %0.2f',f.r));


subplot(2,5,7)
scatter(block2_gamma038_WithNull.dualpairs(:,1),block2_gamma038_WithNull.dualpairs(:,2),'MarkerFaceColor',[0,0.5,1],'MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
hold on;
f = showfit(ezfit(block2_gamma038_WithNull.dualpairs(:,1),block2_gamma038_WithNull.dualpairs(:,2),'affine'),'fitcolor','k','dispeqboxmode','off');
%set(gca, 'XScale','log','YScale','log');
axis square;
text(10e2,10e3,sprintf('r = %0.2f',f.r));

subplot(2,5,8)
scatter(block2_gamma038_WithNull.dualtriple(:,1),block2_gamma038_WithNull.dualtriple(:,2),'MarkerFaceColor',[0,0.5,1],'MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
hold on;
f = showfit(ezfit(block2_gamma038_WithNull.dualtriple(:,1),block2_gamma038_WithNull.dualtriple(:,2),'affine'),'fitcolor','k','dispeqboxmode','off');
%set(gca, 'XScale','log','YScale','log');
axis square;
text(10e2,10e3,sprintf('r = %0.2f',f.r));

subplot(2,5,9)
scatter(block2_gamma038_WithNull.dualquart(:,1),block2_gamma038_WithNull.dualquart(:,2),'MarkerFaceColor',[0,0.5,1],'MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
hold on;
f = showfit(ezfit(block2_gamma038_WithNull.dualquart(:,1),block2_gamma038_WithNull.dualquart(:,2),'affine'),'fitcolor','k','dispeqboxmode','off');
%set(gca, 'XScale','log','YScale','log');
axis square;
text(10e2,10e3,sprintf('r = %0.2f',f.r));

subplot(2,5,10)
scatter(block2_gamma038_WithNull.dualpent(:,1),block2_gamma038_WithNull.dualpent(:,2),'MarkerFaceColor',[0,0.5,1],'MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
hold on;
f = showfit(ezfit(block2_gamma038_WithNull.dualpent(:,1),block2_gamma038_WithNull.dualpent(:,2),'affine'),'fitcolor','k','dispeqboxmode','off');
%set(gca, 'XScale','log','YScale','log');
axis square;
text(10e2,10e3,sprintf('r = %0.2f',f.r));


figure;


subplot(2,5,2)
scatter(log10(block1_gamma038_WithNull.dualpairs(:,1)),log10(block1_gamma038_WithNull.dualpairs(:,2)),'MarkerFaceColor',[1,0.5,0],'MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
hold on;
f = showfit(ezfit(log10(block1_gamma038_WithNull.dualpairs(:,1)),log10(block1_gamma038_WithNull.dualpairs(:,2)),'affine'),'fitcolor','k','dispeqboxmode','off');
%set(gc);
axis square;
text(4,4,sprintf('r = %0.2f',f.r));


subplot(2,5,3)
scatter(log10(block1_gamma038_WithNull.dualtriple(:,1)),log10(block1_gamma038_WithNull.dualtriple(:,2)),'MarkerFaceColor',[1,0.5,0],'MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
hold on;
f = showfit(ezfit(log10(block1_gamma038_WithNull.dualtriple(:,1)),log10(block1_gamma038_WithNull.dualtriple(:,2)),'affine'),'fitcolor','k','dispeqboxmode','off');
%set(gca, 'XScale','log','YScale','log');
axis square;
text(4,4,sprintf('r = %0.2f',f.r));

subplot(2,5,4)
scatter(log10(block1_gamma038_WithNull.dualquart(:,1)),log10(block1_gamma038_WithNull.dualquart(:,2)),'MarkerFaceColor',[1,0.5,0],'MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
hold on;
f = showfit(ezfit(log10(block1_gamma038_WithNull.dualquart(:,1)),log10(block1_gamma038_WithNull.dualquart(:,2)),'affine'),'fitcolor','k','dispeqboxmode','off');
%set(gca, 'XScale','log','YScale','log');
axis square;
text(4,4,sprintf('r = %0.2f',f.r));

subplot(2,5,5)
scatter(log10(block1_gamma038_WithNull.dualpent(:,1)),log10(block1_gamma038_WithNull.dualpent(:,2)),'MarkerFaceColor',[1,0.5,0],'MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
hold on;
f = showfit(ezfit(log10(block1_gamma038_WithNull.dualpent(:,1)),log10(block1_gamma038_WithNull.dualpent(:,2)),'affine'),'fitcolor','k','dispeqboxmode','off');
%set(gca, 'XScale','log','YScale','log');
axis square;
text(4,4,sprintf('r = %0.2f',f.r));



% subplot(2,5,6)
% scatter(block2_gamma038_WithNull.dualpairs(:,1),block2_gamma038_WithNull.dualpairs(:,2),'MarkerFaceColor',[0,0.5,1],'MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
% hold on;
% f = showfit(ezfit(block2_gamma038_WithNull.dualpairs(:,1),block2_gamma038_WithNull.dualpairs(:,2),'affine'),'fitcolor','k','dispeqboxmode','off');
% %set(gca, 'XScale','log','YScale','log');
% axis square;
% text(10e2,10e3,sprintf('r = %0.2f',f.r));


subplot(2,5,7)
scatter(log10(block2_gamma038_WithNull.dualpairs(:,1)),log10(block2_gamma038_WithNull.dualpairs(:,2)),'MarkerFaceColor',[0,0.5,1],'MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
hold on;
f = showfit(ezfit(log10(block2_gamma038_WithNull.dualpairs(:,1)),log10(block2_gamma038_WithNull.dualpairs(:,2)),'affine'),'fitcolor','k','dispeqboxmode','off');
%set(gca, 'XScale','log','YScale','log');
axis square;
text(4,4,sprintf('r = %0.2f',f.r));

subplot(2,5,8)
scatter(log10(block2_gamma038_WithNull.dualtriple(:,1)),log10(block2_gamma038_WithNull.dualtriple(:,2)),'MarkerFaceColor',[0,0.5,1],'MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
hold on;
f = showfit(ezfit(log10(block2_gamma038_WithNull.dualtriple(:,1)),log10(block2_gamma038_WithNull.dualtriple(:,2)),'affine'),'fitcolor','k','dispeqboxmode','off');
%set(gca, 'XScale','log','YScale','log');
axis square;
text(4,4,sprintf('r = %0.2f',f.r));

subplot(2,5,9)
scatter(log10(block2_gamma038_WithNull.dualquart(:,1)),log10(block2_gamma038_WithNull.dualquart(:,2)),'MarkerFaceColor',[0,0.5,1],'MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
hold on;
f = showfit(ezfit(log10(block2_gamma038_WithNull.dualquart(:,1)),log10(block2_gamma038_WithNull.dualquart(:,2)),'affine'),'fitcolor','k','dispeqboxmode','off');
%set(gca, 'XScale','log','YScale','log');
axis square;
text(4,4,sprintf('r = %0.2f',f.r));

subplot(2,5,10)
scatter(log10(block2_gamma038_WithNull.dualpent(:,1)),log10(block2_gamma038_WithNull.dualpent(:,2)),'MarkerFaceColor',[0,0.5,1],'MarkerFaceAlpha',0.2,'MarkerEdgeColor','none');
hold on;
f = showfit(ezfit(log10(block2_gamma038_WithNull.dualpent(:,1)),log10(block2_gamma038_WithNull.dualpent(:,2)),'affine'),'fitcolor','k','dispeqboxmode','off');
%set(gca, 'XScale','log','YScale','log');
axis square;
text(4,4,sprintf('r = %0.2f',f.r));

%%
figure;
subplot(4,4,1)
histogram(block1_gamma038_WithNull.b1tob1size,'binWidth',500,'FaceColor',[1,0.5,0],'EdgeColor','none');
hold on
histogram(block2_gamma038_WithNull.b2tob2size,'binWidth',500,'FaceColor',[0,0.5,1],'EdgeColor','none');
legend({'b1-->b1','b2-->b2'});
xlabel('voxels');
ylabel('count');
axis square;
box off;
offsetAxes(gca);

subplot(4,4,2)
histogram(log10(block1_gamma038_WithNull.b1tob1size),50,'FaceColor',[1,0.5,0],'EdgeColor','none');
hold on
histogram(log10(block2_gamma038_WithNull.b2tob2size),50,'FaceColor',[0,0.5,1],'EdgeColor','none');
legend({'b1-->b1','b2-->b2'});
xlabel('log voxels');
ylabel('count');
axis square;
box off;
set(gca,'XTickLabels',[10^2,10^3,10^4]);
offsetAxes(gca);

subplot(4,4,3)
histogram(block1_gamma038_WithNull.b1tob2size,'binWidth',500,'EdgeColor',[1,0.5,0],'DisplayStyle','stairs','LineWidth',2);
hold on;
histogram(block2_gamma038_WithNull.b2tob1size,'binWidth',500,'EdgeColor',[0,0.5,1],'DisplayStyle','stairs','LineWidth',2);
legend({'b1-->b2','b2-->b1'});
axis square;
box off;
offsetAxes(gca);

subplot(4,4,4)
histogram(log10(block1_gamma038_WithNull.b1tob2size),20,'EdgeColor',[1,0.5,0],'DisplayStyle','stairs','LineWidth',2);
hold on;
histogram(log10(block2_gamma038_WithNull.b2tob1size),20,'EdgeColor',[0,0.5,1],'DisplayStyle','stairs','LineWidth',2);
legend({'b1-->b2','b2-->b1'});
axis square;
xlabel('log voxels');
set(gca,'XTickLabels',[10^2,10^3,10^4]);
box off;
offsetAxes(gca);



subplot(4,4,5)
heatmap([mean(block1_gamma038_WithNull.b1tob1size),mean(block2_gamma038_WithNull.b2tob1size);mean(block1_gamma038_WithNull.b1tob2size),mean(block2_gamma038_WithNull.b2tob2size)],...
    'ColorbarVisible','off','XData',{'b1','b2'},'YData',{'b1','b2'},'FontSize',10);
title('mean synapse size (vx)');

%% size vs numbers

block1_gamma038_WithNull.img_2 = PositionDensity(block1_gamma038_WithNull.cellIDs,2,0);
block1_gamma038_WithNull.img_3 = PositionDensity(block1_gamma038_WithNull.cellIDs,3,0);
block1_gamma038_WithNull.img_4 = PositionDensity(block1_gamma038_WithNull.cellIDs,4,0);
block1_gamma038_WithNull.img_5 = PositionDensity(block1_gamma038_WithNull.cellIDs,5,0);


block2_gamma038_WithNull.img_2 = PositionDensity(block2_gamma038_WithNull.cellIDs,2,0);
block2_gamma038_WithNull.img_3 = PositionDensity(block2_gamma038_WithNull.cellIDs,3,0);
block2_gamma038_WithNull.img_4 = PositionDensity(block2_gamma038_WithNull.cellIDs,4,0);
block2_gamma038_WithNull.img_5 = PositionDensity(block2_gamma038_WithNull.cellIDs,5,0);

% reshaped matrix

block1_gamma038_WithNull.img_2reshape = reshape(cell2mat(block1_gamma038_WithNull.img_2),[],1);
block2_gamma038_WithNull.img_2reshape = reshape(cell2mat(block2_gamma038_WithNull.img_2),[],1);

block1_gamma038_WithNull.img_3reshape = reshape(cell2mat(block1_gamma038_WithNull.img_3),[],1);
block2_gamma038_WithNull.img_3reshape = reshape(cell2mat(block2_gamma038_WithNull.img_3),[],1);

block1_gamma038_WithNull.img_4reshape = reshape(cell2mat(block1_gamma038_WithNull.img_4),[],1);
block2_gamma038_WithNull.img_4reshape = reshape(cell2mat(block2_gamma038_WithNull.img_4),[],1);

block1_gamma038_WithNull.img_5reshape = reshape(cell2mat(block1_gamma038_WithNull.img_5),[],1);
block2_gamma038_WithNull.img_5reshape = reshape(cell2mat(block2_gamma038_WithNull.img_5),[],1);

%
block1_gamma038_WithNull.img_2_90 = sum(block1_gamma038_WithNull.img_2reshape == 90);
block2_gamma038_WithNull.img_2_90 = sum(block2_gamma038_WithNull.img_2reshape == 90);

block1_gamma038_WithNull.img_3_90 = sum(block1_gamma038_WithNull.img_3reshape == 90);
block2_gamma038_WithNull.img_3_90 = sum(block2_gamma038_WithNull.img_3reshape == 90);

block1_gamma038_WithNull.img_4_90 = sum(block1_gamma038_WithNull.img_4reshape == 90);
block2_gamma038_WithNull.img_4_90 = sum(block2_gamma038_WithNull.img_4reshape == 90);

block1_gamma038_WithNull.img_5_90 = sum(block1_gamma038_WithNull.img_5reshape == 90);
block2_gamma038_WithNull.img_5_90 = sum(block2_gamma038_WithNull.img_5reshape == 90);
%%
% subplot(4,4,1)
% 
% plot([1,2,3,4],[block1_gamma038_WithNull.img_2_90/numel(block1_gamma038_WithNull.img_2reshape);...
%                 block1_gamma038_WithNull.img_3_90/numel(block1_gamma038_WithNull.img_3reshape);...
%                 block1_gamma038_WithNull.img_4_90/numel(block1_gamma038_WithNull.img_4reshape);...
%                 block1_gamma038_WithNull.img_5_90/numel(block1_gamma038_WithNull.img_5reshape)],'-ko');
% hold on;
% ylabel('Eyepos pixel fraction (mod1)');
% yyaxis right
% plot([1,2,3,4],[block2_gamma038_WithNull.img_2_90/numel(block2_gamma038_WithNull.img_2reshape);...
%                 block2_gamma038_WithNull.img_3_90/numel(block2_gamma038_WithNull.img_3reshape);...
%                 block2_gamma038_WithNull.img_4_90/numel(block2_gamma038_WithNull.img_4reshape);...
%                 block2_gamma038_WithNull.img_5_90/numel(block2_gamma038_WithNull.img_5reshape)],'-ro');
% set(gca,'YColor','r');
% ylabel('Eyepos pixel fraction (mod2)');
% set(gca,'XTickLabels',{'~6x6x8','~9x9x12','~12x12x16','~16x16x20'},'XTickLabelRotation',45);
% axis square;
% box off;

subplot(4,4,1)
plot([1,2,3,4],[block1_gamma038_WithNull.img_2_90, block1_gamma038_WithNull.img_3_90,...
    block1_gamma038_WithNull.img_4_90,block1_gamma038_WithNull.img_5_90],'-o','color',[1,0.5,0],'LineWidth',2);
hold on
plot([1,2,3,4],[block2_gamma038_WithNull.img_2_90, block2_gamma038_WithNull.img_3_90,...
    block2_gamma038_WithNull.img_4_90,block2_gamma038_WithNull.img_5_90],'-o','color',[0,0.5,1],'LineWidth',2);
set(gca,'XTickLabels',{'~6x6x8','~9x9x12','~12x12x16','~16x16x20'},'XTickLabelRotation',45);
legend({'mod1','mod2'});
ylabel('eyepos neurons');
 axis square;
 box off;



%% positionNess

optimalScale = 2;
% control
block1_gamma038_WithNull.img_2 = PositionDensity(block1_gamma038_WithNull.cellIDs,optimalScale,0);
block2_gamma038_WithNull.img_2 = PositionDensity(block2_gamma038_WithNull.cellIDs,optimalScale,0);

% jitter 5um
block1_gamma038_WithNull.img_2_5 = PositionDensity(block1_gamma038_WithNull.cellIDs,optimalScale,5);
block2_gamma038_WithNull.img_2_5 = PositionDensity(block2_gamma038_WithNull.cellIDs,optimalScale,5);

block1_gamma038_WithNull.img_2_5reshape = reshape(cell2mat(block1_gamma038_WithNull.img_2_5),[],1);
block2_gamma038_WithNull.img_2_5reshape = reshape(cell2mat(block2_gamma038_WithNull.img_2_5),[],1);

% jitter 10um

block1_gamma038_WithNull.img_2_10 = PositionDensity(block1_gamma038_WithNull.cellIDs,optimalScale,10);
block2_gamma038_WithNull.img_2_10 = PositionDensity(block2_gamma038_WithNull.cellIDs,optimalScale,10);

block1_gamma038_WithNull.img_2_10reshape = reshape(cell2mat(block1_gamma038_WithNull.img_2_10),[],1);
block2_gamma038_WithNull.img_2_10reshape = reshape(cell2mat(block2_gamma038_WithNull.img_2_10),[],1);


% jitter 5um

block1_gamma038_WithNull.img_2_15 = PositionDensity(block1_gamma038_WithNull.cellIDs,optimalScale,15);
block2_gamma038_WithNull.img_2_15 = PositionDensity(block2_gamma038_WithNull.cellIDs,optimalScale,15);

block1_gamma038_WithNull.img_2_15reshape = reshape(cell2mat(block1_gamma038_WithNull.img_2_15),[],1);
block2_gamma038_WithNull.img_2_15reshape = reshape(cell2mat(block2_gamma038_WithNull.img_2_15),[],1);
%%

subplot(4,4,2)
plot([1,2,3,4],[sum(block1_gamma038_WithNull.img_2reshape == 90);...
                sum(block1_gamma038_WithNull.img_2_5reshape  == 90);...
                sum(block1_gamma038_WithNull.img_2_10reshape ==90);...
                sum(block1_gamma038_WithNull.img_2_15reshape ==90)],'-o','color',[1,0.5,0],'LineWidth',2);    
ylabel('eyepos neurons')
hold on;
plot([1,2,3,4],[sum(block2_gamma038_WithNull.img_2reshape == 90);...
                sum(block2_gamma038_WithNull.img_2_5reshape  == 90);...
                sum(block2_gamma038_WithNull.img_2_10reshape ==90);...
                sum(block2_gamma038_WithNull.img_2_15reshape ==90)],'-o','color',[0,0.5,1],'LineWidth',2);
set(gca,'XTickLabels',{'centered','5\mu','10\mu','15\mu'},'XTickLabelRotation',35);
box off;
title('@ ~6x6x8');
axis square;

subplot(4,4,3)        
hold on;
plot([1,2,3,4],[sum(block2_gamma038_WithNull.img_2reshape == 90)/sum(block1_gamma038_WithNull.img_2reshape == 90);...
                sum(block2_gamma038_WithNull.img_2_5reshape  == 90)/sum(block1_gamma038_WithNull.img_2_5reshape  == 90);...
                sum(block2_gamma038_WithNull.img_2_10reshape ==90)/sum(block1_gamma038_WithNull.img_2_10reshape ==90);...
                sum(block2_gamma038_WithNull.img_2_15reshape ==90)/ sum(block1_gamma038_WithNull.img_2_15reshape ==90)],'-or');
set(gca,'XTickLabels',{'centered','5\mu','10\mu','15\mu'},'XTickLabelRotation',35);
ylabel('mod2/mod1');
box off;
title('@ ~6x6x8');

axis square;

%%
subplot(4,4,2)

plot([1,2,3,4],[sum(block1_gamma038_WithNull.img_2reshape == 90)/numel(block1_gamma038_WithNull.img_2reshape);...
                sum(block1_gamma038_WithNull.img_2_5reshape  == 90)/numel(block1_gamma038_WithNull.img_2_5reshape);...
                sum(block1_gamma038_WithNull.img_2_10reshape ==90)/ numel(block1_gamma038_WithNull.img_2_10reshape);...
                sum(block1_gamma038_WithNull.img_2_15reshape ==90)/numel(block1_gamma038_WithNull.img_2_15reshape)],'-ok');          
ylabel('eyepos pixel fraction');
yyaxis right;

plot([1,2,3,4],[sum(block2_gamma038_WithNull.img_2reshape == 90)/numel(block2_gamma038_WithNull.img_2reshape);...
                sum(block2_gamma038_WithNull.img_2_5reshape  == 90)/numel(block2_gamma038_WithNull.img_2_5reshape);...
                sum(block2_gamma038_WithNull.img_2_10reshape ==90)/ numel(block2_gamma038_WithNull.img_2_10reshape);...
                sum(block2_gamma038_WithNull.img_2_15reshape ==90)/numel(block2_gamma038_WithNull.img_2_15reshape)],'-or');
set(gca,'YColor','r');
set(gca,'XTickLabels',{'centered','5\mu','10\mu','15\mu'},'XTickLabelRotation',35);
box off;
title('@ ~6x6x8');

subplot(4,4,3)

plot([1,2,3,4],[(sum(block2_gamma038_WithNull.img_2reshape == 90)/numel(block2_gamma038_WithNull.img_2reshape))/...
                (sum(block1_gamma038_WithNull.img_2reshape == 90)/numel(block1_gamma038_WithNull.img_2reshape));...
                 (sum(block2_gamma038_WithNull.img_2_5reshape  == 90)/numel(block2_gamma038_WithNull.img_2_5reshape))/...
                 (sum(block1_gamma038_WithNull.img_2_5reshape  == 90)/numel(block1_gamma038_WithNull.img_2_5reshape));...
                 (sum(block2_gamma038_WithNull.img_2_10reshape ==90)/ numel(block2_gamma038_WithNull.img_2_10reshape))/...
                 (sum(block1_gamma038_WithNull.img_2_10reshape ==90)/ numel(block1_gamma038_WithNull.img_2_10reshape));...
                 (sum(block2_gamma038_WithNull.img_2_15reshape ==90)/numel(block2_gamma038_WithNull.img_2_15reshape))/...
                 (sum(block1_gamma038_WithNull.img_2_15reshape ==90)/numel(block1_gamma038_WithNull.img_2_15reshape))],'-o');
set(gca,'XTickLabels',{'centered','5\mu','10\mu','15\mu'},'XTickLabelRotation',35);
box off;
ylabel('mod2/mod1');

%%
% 
% patchSize = [2*2*2*0.798,2*2*2*0.798,2/2*2*2;...
%              2*3*2*0.798,2*3*2*0.798,3/2*2*2;...
%              2*4*2*0.798,2*4*2*0.798,4/2*2*2;...
%              2*5*2*0.798,2*5*2*0.798,5/2*2*2];

% subplot(4,4,1)
% plot([1,2,3,4],[length(block2_gamma038_WithNull.img_2(block2_gamma038_WithNull.img_2 == 90))/length(block1_gamma038_WithNull.img_2(block1_gamma038_WithNull.img_2 == 90)),...
%    length(block2_gamma038_WithNull.img_3(block2_gamma038_WithNull.img_3 == 90))/length(block1_gamma038_WithNull.img_3(block1_gamma038_WithNull.img_3 == 90)),...
%    length(block2_gamma038_WithNull.img_4(block2_gamma038_WithNull.img_4 == 90))/length(block1_gamma038_WithNull.img_4(block1_gamma038_WithNull.img_4 == 90)),...
%    length(block2_gamma038_WithNull.img_5(block2_gamma038_WithNull.img_5 == 90))/length(block1_gamma038_WithNull.img_5(block1_gamma038_WithNull.img_5 == 90))],'-o');
% set(gca,'XTickLabels',{'~6x6x8','~9x9x12','~12x12x16','~16x16x20'},'XTickLabelRotation',45);
% box off;
% daspect([2,2,1]);
% xlabel('patch size (\mum)')
% ylabel('mod2/mod1 cells')

subplot(4,4,2)
plot([1,2,3,4],...
    [length(block2_gamma038_WithNull.img_2(block2_gamma038_WithNull.img_2 == 90))/length(block1_gamma038_WithNull.img_2(block1_gamma038_WithNull.img_2 == 90)),...
    length(block2_gamma038_WithNull.img_2_5(block2_gamma038_WithNull.img_2_5 == 90))/length(block1_gamma038_WithNull.img_2_5(block1_gamma038_WithNull.img_2_5 == 90)),...
    length(block2_gamma038_WithNull.img_2_10(block2_gamma038_WithNull.img_2_10 == 90))/length(block1_gamma038_WithNull.img_2_10(block1_gamma038_WithNull.img_2_10 == 90)),...
    length(block2_gamma038_WithNull.img_2_15(block2_gamma038_WithNull.img_2_15 == 90))/length(block1_gamma038_WithNull.img_2(block1_gamma038_WithNull.img_2_15 == 90))],...
    '-o');
set(gca,'XTickLabels',{'centered','5\mu','10\mu','15\mu'},'XTickLabelRotation',35);
box off;
daspect([1,2,1]);
xlabel('centroid location')
ylabel('mod2/mod1 cells')
title('@ ~6x6x8')

subplot(4,4,3)

plot([1,2,3,4],...
    [length(block1_gamma038_WithNull.img_2(block1_gamma038_WithNull.img_2 == 90)),...
    length(block1_gamma038_WithNull.img_2_5(block1_gamma038_WithNull.img_2_5 == 90)),...
    length(block1_gamma038_WithNull.img_2_10(block1_gamma038_WithNull.img_2_10 == 90)),...
    length(block1_gamma038_WithNull.img_2(block1_gamma038_WithNull.img_2_15 == 90))],...
    '-o');
hold on;
plot([1,2,3,4],...
    [length(block2_gamma038_WithNull.img_2(block2_gamma038_WithNull.img_2 == 90)),...
    length(block2_gamma038_WithNull.img_2_5(block2_gamma038_WithNull.img_2_5 == 90)),...
    length(block2_gamma038_WithNull.img_2_10(block2_gamma038_WithNull.img_2_10 == 90)),...
    length(block2_gamma038_WithNull.img_2_15(block2_gamma038_WithNull.img_2_15 == 90))],...
        '-o');    
set(gca,'XTickLabels',{'centered','5\mu','10\mu','15\mu'},'XTickLabelRotation',35);
box off;
daspect([1,4,1]);
xlabel('centroid location')
ylabel('number of cells')
legend({'mod1','mod2'})


% block1_gamma038_WithNull_gamma09_WithNull.img = PositionDensity(block1_gamma038_WithNull_gamma09_WithNull.cellIDs);
% block2_gamma038_WithNull_gamma09_WithNull.img = PositionDensity(block2_gamma038_WithNull_gamma09_WithNull.cellIDs);
% block3_gamma09_WithNull.img = PositionDensity(block3_gamma09_WithNull.cellIDs);
% 


figure(1)
for i = 1:length(block1_gamma038_WithNull.cellIDs)
    subplot(13,13,i);
    imagesc(block1_gamma038_WithNull.img_2(:,:,i));
    axis off;
    axis square;
    box on;
    title(block1_gamma038_WithNull.cellIDs(i),'FontSize',8);
end
sgtitle('block1_gamma038_WithNull');
%colormap(colorcet('L19'));

figure(2);
for i = 1:length(block2_gamma038_WithNull.cellIDs)
    subplot(13,13,i);
    imagesc(block2_gamma038_WithNull.img_2(:,:,i));
    axis off;
    axis square;
    box on;
    title(block2_gamma038_WithNull.cellIDs(i),'FontSize',8);
end
sgtitle('block2_gamma038_WithNull');
% 
% figure(3);
% for i = 1:length(block1_gamma038_WithNull_gamma09_WithNull.cellIDs)
%     subplot(13,13,i);
%     imagesc(block1_gamma038_WithNull_gamma09_WithNull.img(:,:,i));
%     axis off;
%     axis square;
%     box on;
%     title(block1_gamma038_WithNull_gamma09_WithNull.cellIDs(i),'FontSize',8);
% end
% sgtitle('block1_gamma038_WithNull_gamma09_WithNull');
% 
% figure(4);
% for i = 1:length(block2_gamma038_WithNull_gamma09_WithNull.cellIDs)
%     subplot(13,13,i);
%     imagesc(block2_gamma038_WithNull_gamma09_WithNull.img(:,:,i));
%     axis off;
%     axis square;
%     box on;
%     title(block2_gamma038_WithNull_gamma09_WithNull.cellIDs(i),'FontSize',8);
% end
% sgtitle('block2_gamma038_WithNull_gamma09_WithNull');
% 
% figure(5);
% for i = 1:length(block3_gamma09_WithNull.cellIDs)
%     subplot(13,13,i);
%     imagesc(block3_gamma09_WithNull.img(:,:,i));
%     axis off;
%     axis square;
%     box on;
%     title(block3_gamma09_WithNull.cellIDs(i),'FontSize',8);
% end
% sgtitle('block3_gamma09_WithNull');

%%

figure;

subplot(1,2,1)
polarhistogram(block1_gamma038_WithNull.img_2(:),25,'Facecolor',[1,0.5,0])
hold on
polarhistogram(block2_gamma038_WithNull.img_2(:),25,'Facecolor',[0,0.5,1])


subplot(1,2,2)
polarhistogram(block1_gamma038_WithNull.img_2_10(:),25,'Facecolor',[1,0.5,0])
hold on
polarhistogram(block2_gamma038_WithNull.img_2_10(:),25,'Facecolor',[0,0.5,1])



%%

figure;

subplot(4,4,2)
histogram(block1_gamma038_WithNull.img_2,'BinWidth',5,'EdgeColor',[1,0.5,0],'DisplayStyle','stairs','LineWidth',2);
hold on
histogram(block2_gamma038_WithNull.img_2,'BinWidth',5,'EdgeColor',[0,0.5,1],'DisplayStyle','stairs','LineWidth',2);
axis square;
box off;
xlabel('grey value');
ylabel('freq');


subplot(4,4,5)
histogram(block1_gamma038_WithNull.img_2_5,'BinWidth',5,'EdgeColor',[1,0.5,0],'DisplayStyle','stairs','LineWidth',2);
hold on
histogram(block2_gamma038_WithNull.img_2_5,'BinWidth',5,'EdgeColor',[0,0.5,1],'DisplayStyle','stairs','LineWidth',2);
%set(gca,'YLim',[0,0.04]);
axis square;
box off; 
title('5\mum jitter');

subplot(4,4,6)
histogram(block1_gamma038_WithNull.img_2_10,'BinWidth',5,'EdgeColor',[1,0.5,0],'DisplayStyle','stairs','LineWidth',2);
hold on
histogram(block2_gamma038_WithNull.img_2_10,'BinWidth',5,'EdgeColor',[0,0.5,1],'DisplayStyle','stairs','LineWidth',2);
%set(gca,'YLim',[0,0.04]);
axis square;
box off;
title('10\mum jitter');


subplot(4,4,7)
histogram(block1_gamma038_WithNull.img_2_15,'BinWidth',5,'EdgeColor',[1,0.5,0],'DisplayStyle','stairs','LineWidth',2);
hold on
histogram(block2_gamma038_WithNull.img_2_15,'BinWidth',5,'EdgeColor',[0,0.5,1],'DisplayStyle','stairs','LineWidth',2);
%set(gca,'YLim',[0,0.04]);
axis square;
box off;
title('15\mum jitter');



%%

errorbar([nanmean(block1_gamma038_WithNull.img{2}(:));nanmean(block1_gamma038_WithNull.img{4}(:));nanmean(block1_gamma038_WithNull.img{6}(:));nanmean(block1_gamma038_WithNull.img{8}(:));nanmean(block1_gamma038_WithNull.img{10}(:));nanmean(block1_gamma038_WithNull.img{12}(:))],...
    [nanstd(block1_gamma038_WithNull.img{2}(:));nanstd(block1_gamma038_WithNull.img{4}(:));nanstd(block1_gamma038_WithNull.img{6}(:));nanstd(block1_gamma038_WithNull.img{8}(:));nanstd(block1_gamma038_WithNull.img{10}(:));nanstd(block1_gamma038_WithNull.img{12}(:))]./[157;157;157;157;157;157])

hold on;
errorbar([nanmean(block2_gamma038_WithNull.img{2}(:));nanmean(block2_gamma038_WithNull.img{4}(:));nanmean(block2_gamma038_WithNull.img{6}(:));nanmean(block2_gamma038_WithNull.img{8}(:));nanmean(block2_gamma038_WithNull.img{10}(:));nanmean(block2_gamma038_WithNull.img{12}(:))],...
    [nanstd(block2_gamma038_WithNull.img{2}(:));nanstd(block2_gamma038_WithNull.img{4}(:));nanstd(block2_gamma038_WithNull.img{6}(:));nanstd(block2_gamma038_WithNull.img{8}(:));nanstd(block2_gamma038_WithNull.img{10}(:));nanstd(block2_gamma038_WithNull.img{12}(:))]./[116;116;116;116;116;116])
