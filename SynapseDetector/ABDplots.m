% density plots

clc
clear;
if ismac
    addpath(genpath('/Users/ashwin/Documents/'));
    df = readtable('/Users/ashwin/Google Drive/Zfish/SynapseDetector/04152019.csv');
else
    addpath(genpath('/usr/people/ashwinv/seungmount/research/Ashwin/Scripts'));
    df = readtable('/usr/people/ashwinv/seungmount/research/Ashwin/SynapseDetector/11252018.csv');
end

startup

load ABDr.mat
load ABDc.mat
load ABDIr.mat
load ABDIc.mat
colorSchemes

%%


ABDrCoords = [];
for i = 1:numel(ABDr)
  % ABDrCoords = [ABDrCoords;[ABDr(i).Tree{1}.Y,ABDr(i).Tree{1}.Z]];
    ABDrCoords = [ABDrCoords;ABDr(i).PreSynCoordsTransformed];
end

ABDIrCoords = [];
for i = 1:numel(ABDIr)
   % ABDIrCoords = [ABDIrCoords;[ABDIr(i).Tree{1}.Y, ABDIr(i).Tree{1}.Z]];
    ABDIrCoords = [ABDIrCoords;ABDIr(i).PreSynCoordsTransformed];
end

ABDcCoords = [];
for i = 1:numel(ABDc)
    if ~isempty(ABDc(i).Tree)
   % ABDcCoords = [ABDcCoords;[ABDc(i).Tree{1}.Y, ABDc(i).Tree{1}.Z]];
    ABDcCoords = [ABDcCoords;ABDc(i).PreSynCoordsTransformed];
    end
end

ABDIcCoords = [];
for i = 1:numel(ABDIc)
  %  ABDIcCoords = [ABDIcCoords;[ABDIc(i).Tree{1}.Y, ABDIc(i).Tree{1}.Z]];
    ABDIcCoords = [ABDIcCoords; ABDIc(i).PreSynCoordsTransformed];

end

subplot(4,4,[1,2])
 h = scatter(ABDrCoords(:,2),ABDrCoords(:,3),'.','MarkerFaceColor',ABDcolor,'MarkerEdgeColor',ABDcolor)
 hold on;
 scatter(ABDIrCoords(:,2),ABDIrCoords(:,3),'.','MarkerFaceColor',ABDicolor,'MarkerEdgeColor',ABDicolor)
 scatter(ABDcCoords(:,2),ABDcCoords(:,3),'.','MarkerFaceColor',ABDcolor,'MarkerEdgeColor',ABDcolor)
 scatter(ABDIcCoords(:,2),ABDIcCoords(:,3),'.','MarkerFaceColor',ABDicolor,'MarkerEdgeColor',ABDicolor)
 set(gca,'Ydir','reverse','XDir','reverse');
 daspect([1,1,1]);
 box off;
 
 %%
 ABDrSaccadicSites = []
 
 for i = 1:numel(ABDr)
     ABDrSaccadicSites = [ABDrSaccadicSites; ABDr(i).PreSynCoordsTransformed(ABDr(i).isSaccadic,:)];
 end
 
 ABDcSaccadicSites = [];
 
  for i = 1:numel(ABDr)
     ABDcSaccadicSites = [ABDcSaccadicSites; ABDc(i).PreSynCoordsTransformed(ABDc(i).isSaccadic,:)];
  end
 
  ABDIrSaccadicSites = []
 
 for i = 1:numel(ABDIr)
     ABDIrSaccadicSites = [ABDIrSaccadicSites; ABDIr(i).PreSynCoordsTransformed(ABDIr(i).isSaccadic,:)];
 end
 
 ABDIcSaccadicSites = [];
 
  for i = 1:numel(ABDIc)
     ABDIcSaccadicSites = [ABDIcSaccadicSites; ABDIc(i).PreSynCoordsTransformed(ABDIc(i).isSaccadic,:)];
  end
  
  
scatter(ABDrCoords(:,2),ABDrCoords(:,3),'.','MarkerFaceColor',ABDcolor,'MarkerEdgeColor',ABDcolor);
hold on;
scatter(ABDcCoords(:,2),ABDcCoords(:,3),'.','MarkerFaceColor',ABDcolor,'MarkerEdgeColor',ABDcolor)

scatter(ABDrSaccadicSites(:,2),ABDrSaccadicSites(:,3),'.','MarkerFaceColor',SaccABDcolor,'MarkerEdgeColor',SaccABDcolor);
scatter(ABDcSaccadicSites(:,2),ABDcSaccadicSites(:,3),'.','MarkerFaceColor',SaccABDcolor,'MarkerEdgeColor',SaccABDcolor)

scatter(ABDIrSaccadicSites(:,2),ABDIrSaccadicSites(:,3),'.','MarkerFaceColor',SaccABDcolor,'MarkerEdgeColor',SaccABDcolor);
scatter(ABDIcSaccadicSites(:,2),ABDIcSaccadicSites(:,3),'.','MarkerFaceColor',SaccABDcolor,'MarkerEdgeColor',SaccABDcolor)

 set(gca,'Ydir','reverse','XDir','reverse');
 daspect([1,1,1]);
 box off;
 
 
 
 %% Fraction reconstructed
 
ABDr_CellIDs = [82140,82145,82143,77648,82146,77302,77301, 82194,77705,77710,77305,77672,77300,77709];
ABDc_CellIDs = [82213,77654,77646,82212,77295,81172,77657,82197,77682,77154, 77652,77658,77628,77292,77688,82195,77296];
ABDIr_CellIDs = [78553,77150,77886,78547,77631,77158,77665,78552,77668,77618,77634];
ABDIc_CellIDs = [78574,79051,77148,77625,77144,77692,77641,79066,77640,77643];


AllABD = [ABDr_CellIDs,ABDc_CellIDs,ABDIr_CellIDs,ABDIc_CellIDs];
Allmotor = [ABDr_CellIDs,ABDc_CellIDs];
 
%%
figure;
subplot(1,2,1)
transform_swc_AV([ABDr_CellIDs,ABDc_CellIDs],ABDcolor,[],true,false);
subplot(1,2,2)
transform_swc_AV([ABDIr_CellIDs,ABDIc_CellIDs],ABDicolor,[],true,false);

figure;
transform_swc_AV([ABDr_CellIDs,ABDc_CellIDs],ABDcolor,[],false,false);
transform_swc_AV([ABDIr_CellIDs,ABDIc_CellIDs],ABDicolor,[],false,false);


%%

 for i = 1:size(AllABD,2)
      [A,B] = SynapticPartners(AllABD(i),1,df);
      tempa = numel(A);
      tempb = numel(A(A<1e5));
      Recon(i,:) = [tempa,tempb/tempa];
      clear A;
      clear B;
      clear tempa;
      clear tempb;
 end
 
  figure;
  subplot(4,4,1)
 histogram(Recon(:,2),'FaceColor','k');
 axis square;
 box off;
 
 
%% Mono, bino fraction per neuron.


for i = 1:size(ABDr,2)
    motorPattern = isMotor(ABDr(i).Inputs(ABDr(i).Inputs<1e5),df);
    OSI = ((motorPattern(:,2)+motorPattern(:,3))-(motorPattern(:,4)+motorPattern(:,5))) ./ ...
        ((motorPattern(:,2)+motorPattern(:,3))+(motorPattern(:,4)+motorPattern(:,5)));
    OSIhistABDr(i,:) = histcounts(OSI,-1:0.1:1);
    %clear OSI;
end



for i = 1:size(ABDc,2)
    if ~isempty(ABDc(i).Tree)
    motorPattern = isMotor(ABDc(i).Inputs(ABDc(i).Inputs<1e5),df);
    OSI = ((motorPattern(:,2)+motorPattern(:,3))-(motorPattern(:,4)+motorPattern(:,5))) ./ ...
        ((motorPattern(:,2)+motorPattern(:,3))+(motorPattern(:,4)+motorPattern(:,5)));
    OSIhistABDc(i,:) = histcounts(OSI,-1:0.1:1);
    else
    OSIhistABDc(i,:)  = repmat(NaN,1,20);
    end
    clear OSI;
end

for i = 1:size(ABDIr,2)
    if ~isempty(ABDIr(i).Tree)
    motorPattern = isMotor(ABDIr(i).Inputs(ABDIr(i).Inputs<1e5),df);
    OSI = ((motorPattern(:,2)+motorPattern(:,3))-(motorPattern(:,4)+motorPattern(:,5))) ./ ...
        ((motorPattern(:,2)+motorPattern(:,3))+(motorPattern(:,4)+motorPattern(:,5)));
    OSIhistABDIr(i,:) = histcounts(OSI,-1:0.1:1);
    else
    OSIhistABDIr(i,:)  = repmat(NaN,1,20);
    end
    clear OSI;
end


for i = 1:size(ABDIc,2)
    if ~isempty(ABDIc(i).Tree)
    motorPattern = isMotor(ABDIc(i).Inputs(ABDIc(i).Inputs<1e5),df);
    OSI = ((motorPattern(:,2)+motorPattern(:,3))-(motorPattern(:,4)+motorPattern(:,5))) ./ ...
        ((motorPattern(:,2)+motorPattern(:,3))+(motorPattern(:,4)+motorPattern(:,5)));
    OSIhistABDIc(i,:) = histcounts(OSI,-1:0.1:1);
    else
    OSIhistABDIc(i,:)  = repmat(NaN,1,20);
    end
    clear OSI;
end
%%

cmap = colorcet('L19','N',20);

subplot(4,4,[5,9])
heatmap(OSIhistABDr,'ColorScaling','scaledrows','Colormap',cmap,'XDisplayLabels', -0.95:0.1:0.95);

subplot(4,4,[6,10])
heatmap(OSIhistABDc,'ColorScaling','scaledrows','Colormap',cmap,'XDisplayLabels', -0.95:0.1:0.95);

subplot(4,4,[7,11])
heatmap(OSIhistABDIr,'ColorScaling','scaledrows','Colormap',cmap,'XDisplayLabels', -0.95:0.1:0.95);

subplot(4,4,[8,12])
heatmap(OSIhistABDIc,'ColorScaling','scaledrows','Colormap',cmap,'XDisplayLabels', -0.95:0.1:0.95);

subplot(4,4,13)
errorbar(mean(OSIhistABDr),std(OSIhistABDr),'-o','color',ABDcolor,'LineWidth',2);
box off;
set(gca,'XTickLabels',[-1,0,1]);
offsetAxes(gca);

subplot(4,4,14)
errorbar(nanmean(OSIhistABDc),nanstd(OSIhistABDc),'-o','color',ABDcolor,'LineWidth',2);
box off;
set(gca,'XTickLabels',[-1,0,1]);
offsetAxes(gca);


subplot(4,4,15)
errorbar(mean(OSIhistABDIr),std(OSIhistABDIr),'-o','color',ABDicolor,'LineWidth',2);
box off;
set(gca,'XTickLabels',[-1,0,1]);
offsetAxes(gca);


subplot(4,4,16)
errorbar(mean(OSIhistABDIc),std(OSIhistABDIc),'-o','color',ABDicolor,'LineWidth',2);
box off;
set(gca,'XTickLabels',[-1,0,1]);
offsetAxes(gca);

figure;

ABDinputs = [vertcat(ABDr.Inputs);vertcat(ABDc.Inputs)];
ABDinputs = ABDinputs(ABDinputs<1e5);
ABDinputs = unique(ABDinputs);
ABDmotorpattern = isMotor(ABDinputs,df);
ABDOSI = (sum(ABDmotorpattern(:,2:3),2)- sum(ABDmotorpattern(:,4:5),2))./...
    (sum(ABDmotorpattern(:,2:3),2)+ sum(ABDmotorpattern(:,4:5),2));

ABDIinputs = [vertcat(ABDIr.Inputs);vertcat(ABDIc.Inputs)];
ABDIinputs = ABDIinputs(ABDIinputs<1e5);
ABDIinputs = unique(ABDIinputs);
ABDImotorpattern = isMotor(ABDIinputs,df);
ABDIOSI = (sum(ABDImotorpattern(:,2:3),2)- sum(ABDImotorpattern(:,4:5),2))./...
    (sum(ABDImotorpattern(:,2:3),2)+ sum(ABDImotorpattern(:,4:5),2));

ABDcomplexInputs = unique([ABDinputs;ABDIinputs]);
ABDcomplexMotorPattern = isMotor(ABDcomplexInputs,df);

figure('Position',[100 100 350 300]);
g = gramm('x',sum(ABDcomplexMotorPattern(:,2:3),2),'y',sum(ABDcomplexMotorPattern(:,4:5),2)); 
g.geom_point('alpha',1);
g.geom_abline();
g.stat_cornerhist('aspect',0.5,'edges',[-120:5:120]);
g.set_color_options('map',[0.1,0.1,0.1]);
%g.set_limit_extra( [0.05,0.05],[0.05,0.05]);
%g.axe_property('XLim',[0,250],'YLim',[0,250]);
g.set_text_options('base_size',15);
g.set_title('Abducens complex')
g.set_names('x','Synapses on M','y','Synapses on I');
g.draw();
g.export('file_name','Abducens','export_path','/Users/ashwin/Desktop','file_type','png')




figure;

subplot(4,4,1)

shadedErrorBar(-0.95:0.1:0.95,nanmean([OSIhistABDr;OSIhistABDc]),nanstd([OSIhistABDr;OSIhistABDc]),'lineprops',{'Color',ABDcolor,'Linewidth',2});
hold on;
shadedErrorBar(-0.95:0.1:0.95,mean([OSIhistABDIr;OSIhistABDIc]),std([OSIhistABDIr;OSIhistABDIc]),'lineprops',{'Color',ABDicolor,'Linewidth',2});
axis square;
box off;
set(gca,'XTickLabels',[-1,0,1]);
offsetAxes(gca);
xlabel('OSI');
ylabel('Count');

subplot(4,4,2)
histogram(ABDOSI,-1:0.1:1,'FaceColor',ABDcolor,'EdgeColor','k');
hold on;
histogram(ABDIOSI,-1:0.1:1,'FaceColor',ABDicolor,'EdgeColor','k')
axis square;
box off;
set(gca,'XTickLabels',[-1,0,1]);
offsetAxes(gca);
xlabel('OSI');
ylabel('Count');

subplot(4,4,4)
histogram(ABDOSI,-1:0.1:1,'EdgeColor',ABDcolor,'DisplayStyle','stairs');
hold on;
histogram(ABDIOSI,-1:0.1:1,'EdgeColor',ABDicolor,'DisplayStyle','stairs');
axis square;
box off;
set(gca,'XTickLabels',[-1,0,1]);
offsetAxes(gca);
xlabel('OSI');
ylabel('Count');

%%
allABDinputs =[vertcat(ABDr.Inputs);vertcat(ABDc.Inputs);vertcat(ABDIr.Inputs);vertcat(ABDIc.Inputs)];
allABDinputs = unique(allABDinputs);
allABDinputs = allABDinputs(allABDinputs<1e5);
allABDmotorPattern = isMotor(allABDinputs,df);
allABDmotorOSI = (sum(allABDmotorPattern(:,2:3),2)- sum(allABDmotorPattern(:,4:5),2))./...
    (sum(allABDmotorPattern(:,2:3),2)+ sum(allABDmotorPattern(:,4:5),2));

figure;
subplot(4,4,1)
histogram(allABDmotorOSI,-1:0.1:1,'EdgeColor','k');
axis square;
box off;
set(gca,'XTickLabels',[-1,0,1]);
offsetAxes(gca);
xlabel('OSI');
ylabel('Count');

colormap(colorcet('D2','N',21));
colorbar;

%% Monocular contributions

ABDmonocular.CellIDs = ABDinputs(find(ABDOSI ==1));
ABDmonocular.Saccadic = ABDmonocular.CellIDs(isSaccade(ABDmonocular.CellIDs));
ABDmonocular.Vestibular =  ABDmonocular.CellIDs(isVestibular(ABDmonocular.CellIDs));
ABDmonocular.Integrator =  ABDmonocular.CellIDs(isIntegrator(ABDmonocular.CellIDs));
ABDmonocular.Contra =  ABDmonocular.CellIDs(isContra(ABDmonocular.CellIDs));
ABDmonocular.Rest = setdiff(ABDmonocular.CellIDs,[ABDmonocular.Saccadic;ABDmonocular.Vestibular;ABDmonocular.Integrator;ABDmonocular.Contra]);
ABDmonocular.Rest = setdiff(ABDmonocular.Rest,[ABDr.cellID,ABDc.cellID]);

ABDimonocular.CellIDs = ABDIinputs(find(ABDIOSI ==-1));
ABDimonocular.Saccadic = ABDimonocular.CellIDs(isSaccade(ABDimonocular.CellIDs));
ABDimonocular.Vestibular =  ABDimonocular.CellIDs(isVestibular(ABDimonocular.CellIDs));
ABDimonocular.Integrator =  ABDimonocular.CellIDs(isIntegrator(ABDimonocular.CellIDs));
ABDimonocular.Contra =  ABDimonocular.CellIDs(isContra(ABDimonocular.CellIDs));
ABDimonocular.Rest = setdiff(ABDimonocular.CellIDs,[ABDimonocular.Saccadic;ABDimonocular.Vestibular;ABDimonocular.Integrator;ABDimonocular.Contra]);
ABDimonocular.Rest = setdiff(ABDimonocular.Rest,[ABDIr.cellID,ABDIc.cellID]);




%% ABD pieplots

colors = colormap(cbrewer('qual','Dark2',5));

for i = 1:size(ABDr,2)
ABDrSaccadicNumel(i) = numel(ABDr(i).Saccadic);
ABDrVestibularNumel(i) = numel(ABDr(i).Vestibular);
ABDrIntegratorNumel(i) = numel(ABDr(i).Integrator);
ABDrContraNumel(i) = numel(ABDr(i).Contra);
ABDrRestNumel(i) = numel(ABDr(i).EverythingElse);
ABDrTotalNumel(i) = numel(ABDr(i).Inputs);
end

for i = 1:size(ABDc,2)
ABDcSaccadicNumel(i) = numel(ABDc(i).Saccadic);
ABDcVestibularNumel(i) = numel(ABDc(i).Vestibular);
ABDcIntegratorNumel(i) = numel(ABDc(i).Integrator);
ABDcContraNumel(i) = numel(ABDc(i).Contra);
ABDcRestNumel(i) = numel(ABDc(i).EverythingElse);
ABDcTotalNumel(i) = numel(ABDc(i).Inputs);
end

subplot(4,4,1)
labels = {'Saccadic','Vestibular', 'Integrator', 'Contra', 'Rest'};

pie([sum(ABDrSaccadicNumel)+sum(ABDcSaccadicNumel),sum(ABDrVestibularNumel)+sum(ABDcVestibularNumel),...
    sum(ABDrIntegratorNumel)+sum(ABDcIntegratorNumel),sum(ABDrContraNumel)+sum(ABDcContraNumel),...
    sum(ABDrRestNumel)+sum(ABDcRestNumel)]);
colormap(cbrewer('qual','Dark2',5));
legend(labels,'Location','best');


subplot(4,4,2)
plot(1:numel(ABDr),ABDrSaccadicNumel./ABDrTotalNumel,'-o', 'MarkerFaceColor', colors(1,:));
hold on
plot(1:numel(ABDr),ABDrVestibularNumel./ABDrTotalNumel,'-o', 'MarkerFaceColor', colors(2,:));
plot(1:numel(ABDr),ABDrIntegratorNumel./ABDrTotalNumel,'-o', 'MarkerFaceColor', colors(3,:));
plot(1:numel(ABDr),ABDrContraNumel./ABDrTotalNumel,'-o', 'MarkerFaceColor', colors(4,:));
plot(1:numel(ABDr),ABDrRestNumel./ABDrTotalNumel,'-o', 'MarkerFaceColor', colors(5,:));
set(gca,'YLim',[0,1]);
box off;

subplot(4,4,3)
plot(1:numel(ABDc),ABDcSaccadicNumel./ABDcTotalNumel,'-o', 'MarkerFaceColor', colors(1,:));
hold on
plot(1:numel(ABDc),ABDcVestibularNumel./ABDcTotalNumel,'-o', 'MarkerFaceColor', colors(2,:));
plot(1:numel(ABDc),ABDcIntegratorNumel./ABDcTotalNumel,'-o', 'MarkerFaceColor', colors(3,:));
plot(1:numel(ABDc),ABDcContraNumel./ABDcTotalNumel,'-o', 'MarkerFaceColor', colors(4,:));
plot(1:numel(ABDc),ABDcRestNumel./ABDcTotalNumel,'-o', 'MarkerFaceColor', colors(5,:));
set(gca,'YLim',[0,1]);
box off;


% Interneuons

for i = 1:size(ABDIr,2)
ABDIrSaccadicNumel(i) = numel(ABDIr(i).Saccadic);
ABDIrVestibularNumel(i) = numel(ABDIr(i).Vestibular);
ABDIrIntegratorNumel(i) = numel(ABDIr(i).Integrator);
ABDIrContraNumel(i) = numel(ABDIr(i).Contra);
ABDIrRestNumel(i) = numel(ABDIr(i).EverythingElse);
ABDIrTotalNumel(i) = numel(ABDIr(i).Inputs);
end

for i = 1:size(ABDIc,2)
ABDIcSaccadicNumel(i) = numel(ABDIc(i).Saccadic);
ABDIcVestibularNumel(i) = numel(ABDIc(i).Vestibular);
ABDIcIntegratorNumel(i) = numel(ABDIc(i).Integrator);
ABDIcContraNumel(i) = numel(ABDIc(i).Contra);
ABDIcRestNumel(i) = numel(ABDIc(i).EverythingElse);
ABDIcTotalNumel(i) = numel(ABDIc(i).Inputs);
end

subplot(4,4,5)

pie([sum(ABDIrSaccadicNumel)+sum(ABDIcSaccadicNumel),sum(ABDIrVestibularNumel)+sum(ABDIcVestibularNumel),...
    sum(ABDIrIntegratorNumel)+sum(ABDIcIntegratorNumel),sum(ABDIrContraNumel)+sum(ABDIcContraNumel),...
    sum(ABDIrRestNumel)+sum(ABDIcRestNumel)]);
colormap(cbrewer('qual','Dark2',5));
legend(labels,'Location','best');

subplot(4,4,6)
plot(1:numel(ABDIr),ABDIrSaccadicNumel./ABDIrTotalNumel,'-o', 'MarkerFaceColor', colors(1,:));
hold on
plot(1:numel(ABDIr),ABDIrVestibularNumel./ABDIrTotalNumel,'-o', 'MarkerFaceColor', colors(2,:));
plot(1:numel(ABDIr),ABDIrIntegratorNumel./ABDIrTotalNumel,'-o', 'MarkerFaceColor', colors(3,:));
plot(1:numel(ABDIr),ABDIrContraNumel./ABDIrTotalNumel,'-o', 'MarkerFaceColor', colors(4,:));
plot(1:numel(ABDIr),ABDIrRestNumel./ABDIrTotalNumel,'-o', 'MarkerFaceColor', colors(5,:));
set(gca,'YLim',[0,1]);
box off;

subplot(4,4,7)
plot(1:numel(ABDIc),ABDIcSaccadicNumel./ABDIcTotalNumel,'-o', 'MarkerFaceColor', colors(1,:));
hold on
plot(1:numel(ABDIc),ABDIcVestibularNumel./ABDIcTotalNumel,'-o', 'MarkerFaceColor', colors(2,:));
plot(1:numel(ABDIc),ABDIcIntegratorNumel./ABDIcTotalNumel,'-o', 'MarkerFaceColor', colors(3,:));
plot(1:numel(ABDIc),ABDIcContraNumel./ABDIcTotalNumel,'-o', 'MarkerFaceColor', colors(4,:));
plot(1:numel(ABDIc),ABDIcRestNumel./ABDIcTotalNumel,'-o', 'MarkerFaceColor', colors(5,:));
set(gca,'YLim',[0,1]);
box off;




%% plot mono vs bistratified neurons

ABDc_monostratified = [82212,77654,82213,77295,77646];
ABDc_Bistratified = [77662,77154,77688,77292,77658,77657,77682,77296,82195,77652,82197,81172,77628,82196];

ABDr_monostratified = [82145,77299,77302,82140,82143,82146];
ABDr_bistratified = [82192,77709,77672,82193,77305,77710,77300,77648,82194,77301,77705];

% figure;
% title('ABDc_mono');
% for i = 1:numel(ABDc_monostratified)
% subplot(4,4,i)
% plotTree(ABDc_monostratified(i))
% axis off;
% title(ABDc_monostratified(i));
% end

% figure;
% title('ABDc_Bi')
% for i = 1:numel(ABDc_Bistratified)
% subplot(4,4,i)
% plotTree(ABDc_Bistratified(i))
% axis off;
% title(ABDc_Bistratified(i));
% end

%  figure;
% title('ABDr_mono');
% for i = 1:numel(ABDr_monostratified)
% subplot(4,4,i)
% plotTree(ABDr_monostratified(i))
% axis off;
% title(ABDr_monostratified(i));
% end


 figure;
title('ABDr_Bi');
for i = 1:numel(ABDr_bistratified)
subplot(4,4,i)
plotTree(ABDr_bistratified(i))
axis off;
title(ABDr_bistratified(i));
end

%%

ABDcInputs = vertcat(ABDc.Inputs);
ABDcInputs = ABDcInputs(ABDcInputs<1e5);
ABDcInputs = unique(ABDcInputs);
ABDcInputsRhombomeres = isRhombomere(ABDcInputs);
clear temp;
temp = ABDcInputsRhombomeres.cellID(logical(ABDcInputsRhombomeres.r4));
temp = temp(~isContra(temp)); 
%transform_swc_AV(temp,SaccABDcolor,[],true,true,'R4ABDcInputs-Vel');
transform_swc_AV(temp,SaccABDcolor,[],true,true);

clear temp;
temp = ABDcInputsRhombomeres.cellID(logical(ABDcInputsRhombomeres.r5));
temp = temp(~isContra(temp)); 
%transform_swc_AV(temp,SaccABDcolor,[],true,true,'R5ABDcInputs-Vel');
transform_swc_AV(temp,SaccABDcolor,[],true,true);


% ABDr

ABDrInputs = vertcat(ABDr.Inputs);
ABDrInputs = ABDrInputs(ABDrInputs<1e5);
ABDrInputs = unique(ABDrInputs);
ABDrInputsRhombomeres = isRhombomere(ABDrInputs);
clear temp;
temp = ABDrInputsRhombomeres.cellID(logical(ABDrInputsRhombomeres.r4));
temp = temp(~isContra(temp)); 
%transform_swc_AV(temp,SaccABDcolor,[],true,true,'R4ABDrInputs-Vel');
transform_swc_AV(temp,SaccABDcolor,[],true,true);

clear temp;
temp = ABDrInputsRhombomeres.cellID(logical(ABDrInputsRhombomeres.r5));
temp = temp(~isContra(temp)); 
%transform_swc_AV(temp,SaccABDcolor,[],true,true,'R5ABDrInputs-Vel');
transform_swc_AV(temp,SaccABDcolor,[],true,true);


% both

ABDinputs = [vertcat(ABDr.Inputs);vertcat(ABDc.Inputs)];
ABDinputs = ABDinputs(ABDinputs<1e5);
ABDinputs = unique(ABDinputs);
ABDinputsRhombomeres = isRhombomere(ABDinputs);




%% OSI analysis for ecverythingelse

% ABD
ABDEverythingElse = [vertcat(ABDr.EverythingElse);vertcat(ABDc.EverythingElse)];
ABDEverythingElse = ABDEverythingElse(ABDEverythingElse<1e5);
ABDEverythingElse = unique(ABDEverythingElse);
ABDEverythingElseIpsi = ABDEverythingElse(~isContra(ABDEverythingElse));
ABDEverythingElseIpsiMotorDist = isMotor(ABDEverythingElseIpsi,df);

%ABDi
ABDiEverythingElse = [vertcat(ABDIr.EverythingElse);vertcat(ABDIc.EverythingElse)];
ABDiEverythingElse = ABDiEverythingElse(ABDiEverythingElse<1e5);
ABDiEverythingElse = unique(ABDiEverythingElse);
ABDiEverythingElseIpsi = ABDiEverythingElse(~isContra(ABDiEverythingElse));
ABDiEverythingElseIpsiMotorDist = isMotor(ABDiEverythingElseIpsi,df);


% histogram

subplot(4,4,1)
histogram(sum(ABDEverythingElseIpsiMotorDist(:,2:3),2),0:10:100);
hold on;
histogram(sum(ABDEverythingElseIpsiMotorDist(:,4:5),2),0:10:100);
title('Presynaptic to ABD');

subplot(4,4,2)
histogram(sum(ABDiEverythingElseIpsiMotorDist(:,2:3),2),0:10:100);
hold on;
histogram(sum(ABDiEverythingElseIpsiMotorDist(:,4:5),2),0:10:100);
title('Presynaptic to ABDi');


% find neurons that have SaccadicPos and Integrators as Posts targets 
ABDEverythingElseSaccadicVel = ABDEverythingElseIpsi(isPostSynapseSaccadic(ABDEverythingElseIpsi,df) & ...
    isPostSynapseIntegrator(ABDEverythingElseIpsi,df) & isPostSynapseIBN(ABDEverythingElseIpsi,df));

% find neurons that have SaccadicPos and Integrators as Posts targets 
ABDiEverythingElseSaccadicVel = ABDiEverythingElseIpsi(isPostSynapseSaccadic(ABDiEverythingElseIpsi,df) & ...
    isPostSynapseIntegrator(ABDiEverythingElseIpsi,df) & isPostSynapseIBN(ABDiEverythingElseIpsi,df));




%remove medial integrators 
lateralVSaccadic = [80163 80167 80177 80179 76688 80204 80206 80210 76682];
ABDEverythingElseSaccadicVel = setdiff(ABDEverythingElseSaccadicVel,lateralVSaccadic);
ABDEverythingElseSaccadicVelMotorDist = isMotor(ABDEverythingElseSaccadicVel,df);
ABDiEverythingElseSaccadicVel = setdiff(ABDiEverythingElseSaccadicVel,lateralVSaccadic);
ABDiEverythingElseSaccadicVelMotorDist = isMotor(ABDiEverythingElseSaccadicVel,df);


ABD_ABDi_EverythginElseSaccadicVel = intersect(ABDEverythingElseSaccadicVel,ABDiEverythingElseSaccadicVel);

% 
subplot(4,4,3)
scatter(sum(ABDEverythingElseSaccadicVelMotorDist(:,2:3),2), sum(ABDEverythingElseSaccadicVelMotorDist(:,4:5),2),25,ABDcolor);
hold on;
scatter(sum(ABDiEverythingElseSaccadicVelMotorDist(:,2:3),2), sum(ABDiEverythingElseSaccadicVelMotorDist(:,4:5),2),25,ABDicolor);

subplot(4,4,4)
ABDEverythingElseSaccadicVelOSI = (sum(ABDEverythingElseSaccadicVelMotorDist(:,2:3),2)-sum(ABDEverythingElseSaccadicVelMotorDist(:,4:5),2)) ./ ...
    (sum(ABDEverythingElseSaccadicVelMotorDist(:,2:3),2)+sum(ABDEverythingElseSaccadicVelMotorDist(:,4:5),2));
histogram(ABDEverythingElseSaccadicVelOSI,-1:0.1:1,'FaceColor',ABDcolor);
hold on;
ABDiEverythingElseSaccadicVelOSI = (sum(ABDiEverythingElseSaccadicVelMotorDist(:,2:3),2)-sum(ABDiEverythingElseSaccadicVelMotorDist(:,4:5),2)) ./ ...
    (sum(ABDiEverythingElseSaccadicVelMotorDist(:,2:3),2)+sum(ABDiEverythingElseSaccadicVelMotorDist(:,4:5),2));
histogram(ABDiEverythingElseSaccadicVelOSI,-1:0.1:1,'FaceColor',ABDicolor);
axis square;


% cleaned up and restricterd to only the medial axons.

 ABDPutativeSaccadic.cellIDs = setdiff(ABDEverythingElseSaccadicVel,ABD_ABDi_EverythginElseSaccadicVel);
 ABDiPutativeSaccadic.cellIDs = setdiff(ABDiEverythingElseSaccadicVel,ABD_ABDi_EverythginElseSaccadicVel);
 ABD_ABDi_PutativeSaccadic.cellIDs = ABD_ABDi_EverythginElseSaccadicVel;
 
 ABDPutativeSaccadic.cellIDs = [81833,81750,79838,81868,78886,80739,78643,80829,81797,78628,79736,81809,81839];
 ABDiPutativeSaccadic.cellIDs = [81759,77338,78634,78680,77846,78628,77372,78643,81420,80286];
 
 
 % plot populations
 
figure;
transform_swc_AV( ABDPutativeSaccadic.cellIDs,ABDputSaccColor,[],true,false);
figure;
transform_swc_AV( ABDiPutativeSaccadic.cellIDs ,ABDiputSaccColor,[],true,false);
figure;
transform_swc_AV(ABD_ABDi_EverythginElseSaccadicVel,[0.5,0.5,0.5],[],true,false);


 
 % Where on ABD and ABDi are the Synapses.

 
 [~,ABDPutativeSaccadic.ABDrpathLength] = getABDgradient(ABDr,ABDPutativeSaccadic.cellIDs,true);
 [~,ABDiPutativeSaccadic.ABDrpathLength] = getABDgradient(ABDr,ABDiPutativeSaccadic.cellIDs,true);
 [~,ABD_ABDi_PutativeSaccadic.ABDrpathLength] = getABDgradient(ABDr,ABD_ABDi_PutativeSaccadic.cellIDs,true);
 
 [~,ABDPutativeSaccadic.ABDcpathLength] = getABDgradient(ABDc,ABDPutativeSaccadic.cellIDs,true);
 [~,ABDiPutativeSaccadic.ABDcpathLength] = getABDgradient(ABDc,ABDiPutativeSaccadic.cellIDs,true);
 [~,ABD_ABDi_PutativeSaccadic.ABDcpathLength] = getABDgradient(ABDc,ABD_ABDi_PutativeSaccadic.cellIDs,true);
 
 [~,ABDPutativeSaccadic.ABDIrpathLength] = getABDgradient(ABDIr,ABDPutativeSaccadic.cellIDs,true);
 [~,ABDiPutativeSaccadic.ABDIrpathLength] = getABDgradient(ABDIr,ABDiPutativeSaccadic.cellIDs,true);
 [~,ABD_ABDi_PutativeSaccadic.ABDIrpathLength] = getABDgradient(ABDIr,ABD_ABDi_PutativeSaccadic.cellIDs,true);
 
 [~,ABDPutativeSaccadic.ABDIcpathLength] = getABDgradient(ABDIc,ABDPutativeSaccadic.cellIDs,true);
 [~,ABDiPutativeSaccadic.ABDIcpathLength] = getABDgradient(ABDIc,ABDiPutativeSaccadic.cellIDs,true);
 [~,ABD_ABDi_PutativeSaccadic.ABDIcpathLength] = getABDgradient(ABDIc,ABD_ABDi_PutativeSaccadic.cellIDs,true);
 
ABDPutativeSaccadic.meanABDgradient = nanmean(vertcat(ABDPutativeSaccadic.ABDrpathLength,ABDPutativeSaccadic.ABDcpathLength));
ABDPutativeSaccadic.stdABDgradient = nanstd(vertcat(ABDPutativeSaccadic.ABDrpathLength,ABDPutativeSaccadic.ABDcpathLength));
ABDPutativeSaccadic.numberOfABDneurons = 29;

ABDiPutativeSaccadic.meanABDigradient = nanmean(vertcat(ABDiPutativeSaccadic.ABDIrpathLength,ABDiPutativeSaccadic.ABDIcpathLength));
ABDiPutativeSaccadic.stdABDigradient = nanstd(vertcat(ABDiPutativeSaccadic.ABDIrpathLength,ABDiPutativeSaccadic.ABDIcpathLength));
ABDiPutativeSaccadic.numberOfABDineurons = 21;

ABD_ABDi_PutativeSaccadic.meanABD_ABDi_ABDgradient = nanmean(vertcat(ABD_ABDi_PutativeSaccadic.ABDrpathLength,ABD_ABDi_PutativeSaccadic.ABDcpathLength));
ABD_ABDi_PutativeSaccadic.stdABD_ABDi_ABDgradient = nanstd(vertcat(ABD_ABDi_PutativeSaccadic.ABDrpathLength,ABD_ABDi_PutativeSaccadic.ABDcpathLength));

ABD_ABDi_PutativeSaccadic.meanABD_ABDi_ABDigradient = nanmean(vertcat(ABD_ABDi_PutativeSaccadic.ABDIrpathLength,ABD_ABDi_PutativeSaccadic.ABDIcpathLength));
ABD_ABDi_PutativeSaccadic.stdABD_ABDi_ABDigradient = nanstd(vertcat(ABD_ABDi_PutativeSaccadic.ABDIrpathLength,ABD_ABDi_PutativeSaccadic.ABDIcpathLength));

 
figure;
subplot(4,4,1)
errorbar(ABDPutativeSaccadic.meanABDgradient,ABDPutativeSaccadic.stdABDgradient./sqrt(ABDPutativeSaccadic.numberOfABDneurons),...
    '-o','color',ABDputSaccColor,'LineWidth',2,'MarkerFaceColor','w');
hold on;
errorbar(ABD_ABDi_PutativeSaccadic.meanABD_ABDi_ABDgradient,ABD_ABDi_PutativeSaccadic.stdABD_ABDi_ABDgradient./sqrt(ABDPutativeSaccadic.numberOfABDneurons),...
    '-o','color',[0.5,0.5,0.5],'LineWidth',2,'MarkerFaceColor','w');
axis square;
set(gca,'XTickLabels',[0,0.5,1],'color',[ABDcolor,0.1],'YLim',[0,1]);
xlabel('Norm. pathlength');
ylabel('Norm. count');
box off;
offsetAxes(gca);

subplot(4,4,2)
errorbar(ABDiPutativeSaccadic.meanABDigradient,ABDiPutativeSaccadic.stdABDigradient./sqrt(ABDiPutativeSaccadic.numberOfABDineurons),...
    '-o','color',ABDiputSaccColor,'LineWidth',2,'MarkerFaceColor','w');
hold on;
errorbar(ABD_ABDi_PutativeSaccadic.meanABD_ABDi_ABDigradient,ABD_ABDi_PutativeSaccadic.stdABD_ABDi_ABDigradient ./sqrt(ABDiPutativeSaccadic.numberOfABDineurons),...
    '-o','color',[0.5,0.5,0.5],'LineWidth',2,'MarkerFaceColor','w');
axis square;
set(gca,'XTickLabels',[0,0.5,1],'color',[ABDicolor,0.1],'YLim',[0,1]);
xlabel('Norm. pathlength');
ylabel('Norm. count');
box off;
offsetAxes(gca);

%%

function  plotTree(cellID);

    filePath  = '/Users/ashwin/Google Drive/Zfish/LowEMtoHighEM/SWC_all/consensus-20180920/swc/';
    oldFilePath = filePath;
    fileName = sprintf('%d.swc',cellID);
    
    if exist(fullfile(filePath,fileName))==0
        filePath  = '/Users/ashwin/Google Drive/Zfish/LowEMtoHighEM/SWC_all/consensus-20190415/swc/';
        oldFilePath = '/Users/ashwin/Google Drive/Zfish/LowEMtoHighEM/SWC_all/consensus-20180920/swc/';
    end
    
    if exist(fullfile(filePath,fileName))
        [tree,~,~] = load_tree(fullfile(filePath,fileName));
        reSampleFactor = 5000;
        tree = resample_tree(tree,reSampleFactor);
        plot_tree(tree,[], [], [], [], '-3l');
        view(0,90);
    end
end

 
 
 










