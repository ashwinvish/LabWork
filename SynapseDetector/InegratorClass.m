% Classes by anatomy
clear;
addpath(genpath('/Users/ashwin/Documents/'));
colorSchemes;

colors = cbrewer('qual','Dark2',10);

ALXSaccadic = colorcet('CBTL1','N',5); %  linear-tritanopic_krjcw_5-98_c46_n256
temp2 = colorcet('CBL2','N',5);  %  linear-protanopic-deuteranopic_kbw_5-98_c40_n256

leadColor = ALXSaccadic(3,:); % dark red, CB safe
lagColor = temp2(3,:); % dark blue, CB safe

PartnerColors = colorcet('CBTD1','N',10);
lightRed = PartnerColors(10,:);
lightBlue = PartnerColors(1,:);


startup

if ismac
    addpath(genpath('/Users/ashwin/Documents/'));
    df = readtable('/Users/ashwin/Google Drive/Zfish/SynapseDetector/04152019.csv');
else
    addpath(genpath('/usr/people/ashwinv/seungmount/research/Ashwin/Scripts'));
    df = readtable('/usr/people/ashwinv/seungmount/research/Ashwin/SynapseDetector/11252018.csv');
end

%% Load Saccadic Neurons

confirmedALX = [76181 76201 76187 76184 76192 76197 ];
putativeALX = [76657 76623 76701 76661 76673 76690 76677 76624 77327 77580 ...
    77341 77344 77868 77872 77349 76634 77592 77578 77607 78629 77581 77591 ...
    77602 77821 77844 77845 77822 78406 78667 79022 79033 78404 78441 78421 ...
    79044 79046 79048 80221 78853 79017 79852 78451 79042 80596 80606 78911 ...
    79746 80271 79720 79976 77586 77369 78633 80750 77142 79060 78453 80885 ...
    81423 81661 81683 81792];

putativeALX = putativeALX(logical(isPostSynapseIntegrator(putativeALX,df))|logical(isPreSynapseIntegrator(putativeALX,df)));
motorOut = isMotor(putativeALX',df);
putativeALX = putativeALX(sum(motorOut(:,2:end),2)>0);
ALX.cellIDs = [confirmedALX,putativeALX];
ALX.cellIDs = ALX.cellIDs(isExistReRoot(ALX.cellIDs));

% remove ALX cells that do not synapse onto M or I neurons
ALX.motorDistribution  = isMotor(ALX.cellIDs',df);
ALX.cellIDs = ALX.cellIDs(sum(ALX.motorDistribution(:,2:end),2)>0);

confirmedDBX = [76199 76200 76182 76183 76185 76186 76189 76191 76188 ];
putativeDBX = [77582 77588 77595 77597 77598 77599 77605 79040 76399 76828 ...
    76829 78435 76826 76874 76289 76542 76832 76836 76838 76877 76883 76884 ...
    78400 78429 76383 76867 77594 76830 76876 78440 78438 78447 79045 78434 ...
    78450 76878 77683 77450 78448 ];

putativeDBX = putativeDBX(logical(isPreSynapseIntegrator(putativeDBX,df)));


DBX.cellIDs = [confirmedDBX,putativeDBX];
DBX.cellIDs = DBX.cellIDs(isExistReRoot(DBX.cellIDs));

confirmedBARHL = [76198 76190 76193 76194 76195 76196 ];
putativeBARHL = [78452 80224 78391];

IntegratorNeurons = [confirmedALX,putativeALX,confirmedDBX,putativeDBX,confirmedBARHL,putativeBARHL];
IntegratorNeurons = IntegratorNeurons(isExistReRoot(IntegratorNeurons)); % remove those integrators for which we do not have meshes
%%
for i = 1:numel(IntegratorNeurons)
    if ~isExistReRoot(IntegratorNeurons(i))==0
        %disp(i)
        Int(i) = InputsByClass(IntegratorNeurons(i),df,1);
        % locate and remove self touches among the integrators
        ind = find(ismember(Int(i).Integrator,IntegratorNeurons(i)));
        Int(i).Integrator(ind) = [];
        clear ind;
    end
end

%% ALX OSI


SaccadicProjectingToABDexclusively.cellID;
SaccadicProjectingToABDiexclusively.cellID


%%

ConfirmedNeuorons = find(ismember(IntegratorNeurons,[confirmedALX,confirmedDBX,confirmedBARHL]));

for i =1:numel(ConfirmedNeuorons)
    Conf.cellID(i)  = ConfirmedNeuorons(i);
    Conf.SaccadicNumbers(i) = length(Int(ConfirmedNeuorons(i)).Saccadic);
    Conf.VestibularNumbers(i) = length(Int(ConfirmedNeuorons(i)).Vestibular);
    Conf.IntegratorNumbers(i) = length(Int(ConfirmedNeuorons(i)).Integrator);
    Conf.ContraNumbers(i) = length(Int(ConfirmedNeuorons(i)).Contra);
    Conf.RestNumbers(i) = length(Int(ConfirmedNeuorons(i)).EverythingElse);
end

for i =1:numel(IntegratorNeurons)
    All.cellID(i) = IntegratorNeurons(i);
    All.SaccadicNumbers(i) = length(Int(i).Saccadic);
    All.VestibularNumbers(i) = length(Int(i).Vestibular);
    All.IntegratorNumbers(i) = length(Int(i).Integrator);
end

for i = 1:numel(ALX.cellIDs)
    indexi = find(ismember([Int.cellID]',ALX.cellIDs(i)));
    if ~isempty(indexi)
        %ALX(i) = Int(indexi);
        ALX.SaccadicNumbers(i) = length(Int(indexi).Saccadic);
        ALX.VestibularNumbers(i) = length(Int(indexi).Vestibular);
        ALX.IntegratorNumbers(i) = length(Int(indexi).Integrator);
        ALX.ContraNumbers(i) = length(Int(indexi).Contra);
        ALX.RestNumbers(i) = length(Int(indexi).EverythingElse);
        if isempty(Int(indexi).Origin)
            ALX.Origin(i,:) = [NaN,NaN,NaN];
        else
            ALX.Origin(i,:) = Int(indexi).Origin;
        end
    end
end

ALX.TotalInspiNumber = ALX.SaccadicNumbers+ALX.VestibularNumbers+ALX.IntegratorNumbers+ALX.RestNumbers;
ALX.TotalNumber = ALX.TotalInspiNumber + ALX.ContraNumbers;


for i = 1:numel(DBX.cellIDs)
    indexi = find(ismember([Int.cellID]',DBX.cellIDs(i)));
    if ~isempty(indexi)
       % DBX(i) = Int(indexi);
        DBX.SaccadicNumbers(i) = length(Int(indexi).Saccadic);
        DBX.VestibularNumbers(i) = length(Int(indexi).Vestibular);
        DBX.IntegratorNumbers(i) = length(Int(indexi).Integrator);
        DBX.ContraNumbers(i) = length(Int(indexi).Contra);
        DBX.RestNumbers(i) = length(Int(indexi).EverythingElse);
        DBX.Origin(i,:) = Int(indexi).Origin
    end
end

DBX.TotalInspiNumber = DBX.SaccadicNumbers+DBX.VestibularNumbers+DBX.IntegratorNumbers+DBX.RestNumbers;
DBX.TotalNumber = DBX.TotalInspiNumber + DBX.ContraNumbers;

%% Gradients and cumulative plots

load ABDPutativeSaccadic.mat
load ABDiPutativeSaccadic.mat
load SaccadicProjectingToABDexclusively.mat
load SaccadicProjectingToABDiexclusively.mat
lateralVSaccadic = [80163 80167 80177 80179 76688 80204 80206 80210 76682];

allIntegratorCellIDs = [SaccadicProjectingToABDexclusively.cellID';SaccadicProjectingToABDiexclusively.cellID';lateralVSaccadic';ALX.cellIDs'];


[Integ.VestibularPathLength, Integ.VestibularGradient] = getABDgradient(Int, unique(vertcat(Int.Vestibular)),false,false);
[Integ.SaccadicPathLength, Integ.SaccadicGradient] = getABDgradient(Int,[ABDPutativeSaccadic.cellIDs,ABDiPutativeSaccadic.cellIDs] ,false,false);
[Integ.IntegratorPathLength, Integ.IntegratorGradient] = getABDgradient(Int,allIntegratorCellIDs ,false,false);
[Integ.r456IntegratorPathLength, Integ.r456IntegratorGradient] = getABDgradient(Int,[SaccadicProjectingToABDexclusively.cellID';SaccadicProjectingToABDiexclusively.cellID'],false,false);
[Integ.latIntegratorPathLength, Integ.latIntegratorGradient] = getABDgradient(Int,lateralVSaccadic' ,false,false);
[Integ.r78IntegratorPathLength, Integ.r78IntegratorGradient] = getABDgradient(Int,ALX.cellIDs' ,false,false);


[ALX_RCsorted,ALX.RCsort] = sort(ALX.Origin(:,2));
[DBX_RCsorted,DBX.RCsort] = sort(DBX.Origin(:,2));

for i = 1:numel(ALX.cellIDs)
    indexi = find(ismember([Int.cellID]',ALX.cellIDs(i)));
    if ~isempty(indexi)
        ALX.VestibularPathLenght{i} = Integ.VestibularPathLength{indexi};
        ALX.VestibularGradient(i,:) = Integ.VestibularGradient(indexi,:);
        ALX.SaccadicPathLength{i} = Integ.SaccadicPathLength{indexi};
        ALX.SaccadicGradient(i,:) = Integ.SaccadicGradient(indexi,:);
        ALX.IntegratorPathLength{i} = Integ.IntegratorPathLength{indexi};
        ALX.IntegratorGradient(i,:) = Integ.IntegratorGradient(indexi,:);
        ALX.r456IntegratorPathLength{i} = Integ.r456IntegratorPathLength{indexi};
        ALX.r456IntegratorGradient(i,:) = Integ.r456IntegratorGradient(indexi,:);
        ALX.latIntegratorPathLength{i} = Integ.latIntegratorPathLength{indexi};
        ALX.latIntegratorGradient(i,:) = Integ.latIntegratorGradient(indexi,:);
        ALX.r78IntegratorPathLength{i} = Integ.r78IntegratorPathLength{indexi};
        ALX.r78IntegratorGradient(i,:) = Integ.r78IntegratorGradient(indexi,:);
    end
end

for i = 1:numel(DBX.cellIDs)
    indexi = find(ismember([Int.cellID]',DBX.cellIDs(i)));
    if ~isempty(indexi)
        DBX.VestibularPathLenght{i} = Integ.VestibularPathLength{indexi};
        DBX.VestibularGradient(i,:) = Integ.VestibularGradient(indexi,:);
        DBX.SaccadicPathLength{i} = Integ.SaccadicPathLength{indexi};
        DBX.SaccadicGradient(i,:) = Integ.SaccadicGradient(indexi,:);
        DBX.IntegratorPathLength{i} = Integ.IntegratorPathLength{indexi};
        DBX.IntegratorGradient(i,:) = Integ.IntegratorGradient(indexi,:);
        DBX.r456IntegratorPathLength{i} = Integ.r456IntegratorPathLength{indexi};
        DBX.r456IntegratorGradient(i,:) = Integ.r456IntegratorGradient(indexi,:);
        DBX.latIntegratorPathLength{i} = Integ.latIntegratorPathLength{indexi};
        DBX.latIntegratorGradient(i,:) = Integ.latIntegratorGradient(indexi,:);
        DBX.r78IntegratorPathLength{i} = Integ.r78IntegratorPathLength{indexi};
        DBX.r78IntegratorGradient(i,:) = Integ.r78IntegratorGradient(indexi,:);
    end
end

figure;
subplot(4,4,1)
histogram(cell2mat(ALX.SaccadicPathLength),20,'Normalization','cdf','EdgeColor','k','LineWidth',2,'DisplayStyle','stairs');
hold on;
histogram(cell2mat(ALX.VestibularPathLenght),20,'Normalization','cdf','EdgeColor','r','LineWidth',2,'DisplayStyle','stairs');
histogram(cell2mat(ALX.IntegratorPathLength),20,'Normalization','cdf','EdgeColor','b','LineWidth',2,'DisplayStyle','stairs');
axis square;
title('r78ipsi');
box off;
offsetAxes(gca);

subplot(4,4,2)
histogram(cell2mat(DBX.SaccadicPathLength),20,'Normalization','cdf','EdgeColor','k','LineWidth',2,'DisplayStyle','stairs');
hold on;
histogram(cell2mat(DBX.VestibularPathLenght),20,'Normalization','cdf','EdgeColor','r','LineWidth',2,'DisplayStyle','stairs');
histogram(cell2mat(DBX.IntegratorPathLength),20,'Normalization','cdf','EdgeColor','b','LineWidth',2,'DisplayStyle','stairs');
axis square;
title('r78contra');
box off;
offsetAxes(gca);


figure;
sgtitle('ALX');

saccadicColorMap = ['#ffffff'; '#f5f5f5'; '#eaeaea'; '#e0e0e0'; '#d6d6d6'; '#cccccc'; '#c2c2c2'; '#b8b8b8'; '#aeaeae'; '#a4a4a4'; '#9b9b9b'; '#919191'; '#888888'; '#7e7e7e'; '#757575'; '#6c6c6c'; '#636363'; '#5b5b5b'; '#525252'; '#4a4a4a'];
%saccadicColorMap = colorcet('L1','N',20,'reverse',1);
vestibuarColorMap = ['#ffffff'; '#fff9f4'; '#fff3e9'; '#ffedde'; '#ffe7d3'; '#ffe1c7'; '#ffdabc'; '#ffd4b1'; '#ffcea6'; '#ffc79a'; '#ffc18e'; '#ffba83'; '#ffb477'; '#ffad6a'; '#ffa65e'; '#ff9e50'; '#ff9742'; '#ff8f33'; '#ff8720'; '#ff7f00'];
integratorColorMap = ['#ffffff'; '#f5f8fb'; '#ecf1f8'; '#e2e9f4'; '#d9e2f0'; '#cfdbec'; '#c5d4e9'; '#bccde5'; '#b2c6e1'; '#a8c0dd'; '#9eb9da'; '#94b2d6'; '#8aabd2'; '#80a5ce'; '#769ecb'; '#6b98c7'; '#5f91c3'; '#538bbf'; '#4684bc'; '#377eb8'];


subplot(1,5,1)
heatmap(ALX.SaccadicGradient(ALX.RCsort,:),'Colormap',hex2rgb(saccadicColorMap),'MissingDataColor','w');
title('Sacc');

subplot(1,5,2)
heatmap(ALX.VestibularGradient(ALX.RCsort,:),'Colormap',hex2rgb(vestibuarColorMap),'MissingDataColor','w');
title('Vest');

subplot(1,5,3)
heatmap(ALX.r78IntegratorGradient(ALX.RCsort,:),'Colormap',hex2rgb(integratorColorMap),'MissingDataColor','w');
title('r78');

subplot(1,5,4)
heatmap(ALX.r456IntegratorGradient(ALX.RCsort,:),'Colormap',hex2rgb(integratorColorMap),'MissingDataColor','w');
title('r456');

subplot(1,5,5)
heatmap(ALX.latIntegratorGradient(ALX.RCsort,:),'Colormap',hex2rgb(integratorColorMap),'MissingDataColor','w');
title('r56');


figure;
sgtitle('DBX')

subplot(1,5,1)
heatmap(DBX.SaccadicGradient(DBX.RCsort,:),'Colormap',hex2rgb(saccadicColorMap),'MissingDataColor','w');
title('Sacc');

subplot(1,5,2)
heatmap(DBX.VestibularGradient(DBX.RCsort,:),'Colormap',hex2rgb(vestibuarColorMap),'MissingDataColor','w');
title('Vest');

subplot(1,5,3)
heatmap(DBX.r78IntegratorGradient(DBX.RCsort,:),'Colormap',hex2rgb(integratorColorMap),'MissingDataColor','w');
title('r78');

subplot(1,5,4)
heatmap(DBX.r456IntegratorGradient(DBX.RCsort,:),'Colormap',hex2rgb(integratorColorMap),'MissingDataColor','w');
title('r456');

subplot(1,5,5)
heatmap(DBX.latIntegratorGradient(DBX.RCsort,:),'Colormap',hex2rgb(integratorColorMap),'MissingDataColor','w');
title('r56');

%% Total ALX contribution

ALXtoALX = sum(ALX.r78IntegratorGradient,2)./ALX.TotalNumber';
ALXtoALX(find(ALXtoALX == 16)) = [];

AllIntegratortoALX =(sum(ALX.r78IntegratorGradient,2)+ sum(ALX.r456IntegratorGradient,2)) ./ ALX.TotalNumber' ; 
AllIntegratortoALX (find(AllIntegratortoALX  == 16)) = [];

%%

clear s
clear v
clear saccadic_weights
clear vestibular_weights
clear scaling_proximity

scaling_proximity = ones(40,1)*[10:-1:1];
saccadic_weights = sum(ALX.SaccadicGradient(ALX.RCsort,:).*scaling_proximity,2);
vestibular_weights = sum(ALX.VestibularGradient(ALX.RCsort,:).*scaling_proximity,2);
ALX.pos_connection_index = find(saccadic_weights+vestibular_weights > 0);
s = saccadic_weights(ALX.pos_connection_index);v = vestibular_weights(ALX.pos_connection_index);
ALX.gradient_metric = (v-s)./(v+s);

%ALX.smoothed_gradient_metric = smooth(gradient_metric, 2, 1);

clear s
clear v
clear saccadic_weights
clear vestibular_weights
clear scaling_proximity

scaling_proximity = ones(31,1)*[10:-1:1];
saccadic_weights = sum(DBX.SaccadicGradient(DBX.RCsort,:).*scaling_proximity,2);
vestibular_weights = sum(DBX.VestibularGradient(DBX.RCsort,:).*scaling_proximity,2);
DBX.pos_connection_index = find(saccadic_weights+vestibular_weights > 0);
s = saccadic_weights(DBX.pos_connection_index);v = vestibular_weights(DBX.pos_connection_index);

DBX.gradient_metric = (v-s)./(v+s);

figure;
subplot(4,4,1)
scatter(ALX.Origin(ALX.RCsort(ALX.pos_connection_index),2),smooth(ALX.gradient_metric,5),'o',...
    'MarkerFaceColor',ALXcolor,'MarkerEdgeColor','k','MarkerFaceAlpha',0.4);
hold on;
plot(ALX.Origin(ALX.RCsort(ALX.pos_connection_index),2),smooth(ALX.gradient_metric,5),'LineWidth',2,'Color',ALXcolor);
scatter(DBX.Origin(DBX.RCsort(DBX.pos_connection_index),2),smooth(DBX.gradient_metric, 5),'o',...
    'MarkerFaceColor',DBXcolor,'MarkerEdgeColor','k','MarkerFaceAlpha',0.4)
plot(DBX.Origin(DBX.RCsort(DBX.pos_connection_index),2),smooth(DBX.gradient_metric, 20),'LineWidth',2,'Color',DBXcolor);
legend({'ALX','DBX'});
set(gca,'XLim',[650,800]);
box off;
ylabel('V-S/V+S');
xlabel('RC positon');
view([90,90]);
daspect([40,1,1]);
offsetAxes(gca)

subplot(4,4,2)
ALX.IntegratorGradient  = ALX.r78IntegratorGradient+ALX.r456IntegratorGradient+ALX.latIntegratorGradient;
plot(ALX.Origin(:,2),sum(ALX.IntegratorGradient,2)./ALX.TotalInspiNumber','o',...
    'MarkerFaceColor',ALXcolor,'MarkerEdgeColor','k');
hold on;
plot(ALX.Origin(:,2),sum(ALX.IntegratorGradient,2)./ALX.TotalInspiNumber',...
   'LineWidth',2);
hold on;
box off;
ylabel('Input fraction');
xlabel('RC positon');
set(gca,'YLim',[0,0.5]);
view([90,90]);
%daspect([50,1,1]);
offsetAxes(gca)



%% plots
subplot(4,7,1)
violinPlot(ALX.SaccadicNumbers', 'histOri', 'left',  'showMM', 0, 'color', colors(1,:),'widthDiv',[2,1] ,'showMM',2);
hold on;
violinPlot(DBX.SaccadicNumbers', 'histOri', 'right',  'showMM', 0, 'color', colors(2,:),'widthDiv',[2,2],'showMM',2);
ylabel('Synapses');
set(gca, 'xtick', [0.6 1.4], 'xticklabel', {'ALX', 'DBX'},'ylim',[0,60]);
xlabel('Saccadic axons');

subplot(4,7,3)
violinPlot(ALX.VestibularNumbers', 'histOri', 'left',  'showMM', 0, 'color', colors(1,:),'widthDiv',[2,1],'showMM',2);
hold on;
violinPlot(DBX.VestibularNumbers', 'histOri', 'right',  'showMM', 0, 'color', colors(2,:),'widthDiv',[2,2],'showMM',2);
xlabel('Vestibular axons');
set(gca, 'xtick', [0.6 1.4], 'xticklabel', {'ALX', 'DBX'},'ylim',[0,60]);


subplot(4,7,5)
violinPlot(ALX.IntegratorNumbers', 'histOri', 'left',  'showMM', 0, 'color', colors(1,:),'widthDiv',[2,1],'showMM',2);
hold on;
violinPlot(DBX.IntegratorNumbers', 'histOri', 'right',  'showMM', 0, 'color', colors(2,:),'widthDiv',[2,2],'showMM',2);
xlabel('Integrator axons');
set(gca, 'xtick', [0.6 1.4], 'xticklabel', {'ALX', 'DBX'},'ylim',[0,60]);



subplot(4,7,8)
violinPlot(ALX.SaccadicNumbers', 'histOri', 'left',  'showMM', 0, 'color', colors(1,:),'widthDiv',[2,1] ,'showMM',2,'histOpt',0);
hold on;
violinPlot(DBX.SaccadicNumbers', 'histOri', 'right',  'showMM', 0, 'color', colors(2,:),'widthDiv',[2,2],'showMM',2,'histOpt',0);
ylabel('Synapses');
set(gca, 'xtick', [0.6 1.4], 'xticklabel', {'ALX', 'DBX'});

subplot(4,7,10)
violinPlot(ALX.VestibularNumbers', 'histOri', 'left',  'showMM', 0, 'color', colors(1,:),'widthDiv',[2,1],'showMM',2,'histOpt',0);
hold on;
violinPlot(DBX.VestibularNumbers', 'histOri', 'right',  'showMM', 0, 'color', colors(2,:),'widthDiv',[2,2],'showMM',2,'histOpt',0);


subplot(4,7,12)
violinPlot(ALX.IntegratorNumbers', 'histOri', 'left',  'showMM', 0, 'color', colors(1,:),'widthDiv',[2,1],'showMM',2,'histOpt',0);
hold on;
violinPlot(DBX.IntegratorNumbers', 'histOri', 'right',  'showMM', 0, 'color', colors(2,:),'widthDiv',[2,2],'showMM',2,'histOpt',0);


subplot(4,7,15)
scatter(ALX.Origin(:,2),ALX.SaccadicNumbers,25,colors(1,:));
hold on;
scatter(DBX.Origin(:,2),DBX.SaccadicNumbers,25,colors(2,:));
xlabel('R<-->C');
ylabel('Saccadic Synapses');
daspect([2,1,1]);

subplot(4,7,17)
scatter(ALX.Origin(:,2),ALX.VestibularNumbers,25,colors(1,:));
hold on;
scatter(DBX.Origin(:,2),DBX.VestibularNumbers,25,colors(2,:));
xlabel('R<-->C');
ylabel('Vestibular Synapses');
daspect([2,1,1]);

subplot(4,7,19)
scatter(ALX.Origin(:,2),ALX.IntegratorNumbers,25,colors(1,:));
hold on;
scatter(DBX.Origin(:,2),DBX.IntegratorNumbers,25,colors(2,:));
xlabel('R<-->C');
ylabel('Integrator Synapses');
daspect([4,1,1]);

subplot(4,7,22)
histogram(ALX.Origin(:,2),'FaceColor',colors(1,:));
hold on;
histogram(DBX.Origin(:,2),'FaceColor',colors(2,:));
box off;
%%


%% Pie plots of input/putput fractions


ALXPie = [sum(ALX.SaccadicNumbers),sum(ALX.VestibularNumbers),sum(ALX.IntegratorNumbers),...
    sum(ALX.ContraNumbers),sum(ALX.RestNumbers)];

DBXPie = [sum(DBX.SaccadicNumbers),sum(DBX.VestibularNumbers),sum(DBX.IntegratorNumbers),...
    sum(DBX.ContraNumbers),sum(DBX.RestNumbers)];

ConfPie = [sum(Conf.SaccadicNumbers),sum(Conf.VestibularNumbers),sum(Conf.IntegratorNumbers),...
    sum(Conf.ContraNumbers),sum(Conf.RestNumbers)];

figure;

subplot(3,5,1)
p1 = pie(ALXPie);
txt = {'Saccadic';'Vestibular';'Integrator'; 'Contra';'Rest'};
legend(txt,'Location','bestoutside');
colormap(colors(1:5,:));
title('Ipsi integrators');

subplot(3,5,4)
p2 = pie(DBXPie);
txt = {'Saccadic';'Vestibular';'Integrator'; 'Contra';'Rest'};
colormap(colors(1:5,:));
title('Contra Integrators');

subplot(3,5,6)
p3 = pie(ConfPie);
txt = {'Saccadic';'Vestibular';'Integrator'; 'Contra';'Rest'};
colormap(colors(1:5,:));
title('Confirmed Integrators');

%% Ocular Selectivity

ALX.motorDistribution = isMotor(ALX.cellIDs',df);
ALX.ABDcounts = ALX.motorDistribution(:,2)+ALX.motorDistribution(:,3);
ALX.ABDicounts = ALX.motorDistribution(:,4)+ALX.motorDistribution(:,5);
ALX.NormalizedCounts = (ALX.ABDcounts-ALX.ABDicounts) ./ (ALX.ABDcounts+ALX.ABDicounts);


subplot(4,4,1);
histogram(ALX.NormalizedCounts,-1:0.1:1,'FaceColor',ALXcolor);
hold on;
line([0.5,0.5],[0,20],'color','k','LineStyle',':');
line([-0.5,-0.5],[0,20],'color','k','LineStyle',':');
xlabel('OSI');
axis square;
box off;
offsetAxes;

subplot(4,4,2)
scatter(ALX.ABDcounts,ALX.ABDicounts,'o','MarkerFaceColor',ALXcolor,'MarkerEdgeColor','k');
axis square;
xlabel('Synapses on M');
ylabel('Synapses on I');
offsetAxes(gca);




ALXwithZeroValues = find(ALX.ABDcounts == 0 & ALX.ABDicounts ==0 );

subplot(4,4,3)
histogram2(ALX.ABDcounts,ALX.ABDicounts,...
    0:2:100,0:2:100,'DisplayStyle','tile','ShowEmptyBins','off');
title('ALXs');
box off;
axis square;
colormap(colorcet('L17','reverse',1));

subplot(4,4,4)
histogram2(SaccadicMotorDistribution.ABD(SaccadicProjectingToABD_ABDi_Index),...
    SaccadicMotorDistribution.ABDi(SaccadicProjectingToABD_ABDi_Index),...
    0:5:100,0:5:100,'DisplayStyle','tile','ShowEmptyBins','off');
title('Saccadic_ABD_ABDi')



figure('Position',[100 100 350 300]);
g = gramm('x',ALX.ABDcounts(setdiff(1:66,ALXwithZeroValues)),'y',ALX.ABDicounts(setdiff(1:66,ALXwithZeroValues)));
g.geom_point();
g.geom_abline();
g.stat_cornerhist('location',80,'edges',-100:10:100,'aspect',0.5);
g.set_color_options('map',[0.5,0.5,0.5]);
%g.set_limit_extra( [0.05,0.05],[0.05,0.05]);
%g.axe_property('XLim',[0,250],'YLim',[0,250]);
g.set_text_options('base_size',15);
g.set_names('x','Synapses on M','y','Synapses on I');
g.set_title('Integrators')
g.draw();
g.export('file_name','Integrator','export_path','/Users/ashwin/Desktop','file_type','png')



ALXprojectingToABDIndex = find(ALX.ABDcounts>5 & ALX.ABDicounts<5);
ALXprojectingToABD = ALX.cellIDs(ALXprojectingToABDIndex);

ALXprojectingToABDiIndex = find(ALX.ABDcounts<5 & ALX.ABDicounts>5);
ALXprojectingToABDi = ALX.cellIDs(ALXprojectingToABDiIndex);

ALXprojectingToABD_ABDi_Index = find(ALX.ABDcounts>5 & ALX.ABDicounts>5);
ALXprojectingToABD_ABDi = ALX.cellIDs(ALXprojectingToABD_ABDi_Index);




centerALX = numel(find(ALX.NormalizedCounts< 0.5 & ALX.NormalizedCounts> -0.5 ));
rightALX = numel(find(ALX.NormalizedCounts > 0.5));
leftALX = numel(find(ALX.NormalizedCounts < -0.5));

% plot the ALX populations

rightALXcellIDS = ALX.cellIDs(find(ALX.NormalizedCounts > 0.5));
leftALXcellIDS = ALX.cellIDs(find(ALX.NormalizedCounts < -0.5));

figure;
transform_swc_AV(rightALXcellIDS, ALXcolor,[],true,false);
figure;
transform_swc_AV(leftALXcellIDS, ALXcolor,[],true,false);

%% location along abucnes

load ABDr.mat
load ABDc.mat
load ABDIr.mat
load ABDIc.mat

for i = 1:numel(ABDr)
    temp1 = ismember(ABDr(i).Inputs,ALX.cellIDs');
    if sum(temp1) ~=0 ;
        ALX.ABDrpathLength(i,:) = histcounts(ABDr(i).PathLength(temp1)'./(max(Pvec_tree(ABDr(i).Tree{1}))),0:0.1:1,'Normalization','probability');
    else
        ALX.ABDrpathLength(i,:) = repmat(NaN,1,10);
    end
end

for i = 1:numel(ABDc)
    if ~isempty(ABDc(i).Tree)
        temp1 = ismember(ABDc(i).Inputs,ALX.cellIDs');
        if sum(temp1)~=0 ;
            ALX.ABDcpathLength(i,:) = histcounts(ABDc(i).PathLength(temp1)'./(max(Pvec_tree(ABDc(i).Tree{1}))),0:0.1:1,'Normalization','probability');
        else
            ALX.ABDcpathLength(i,:) = repmat(NaN,1,10);
        end
    end
end

for i = 1:numel(ABDIr)
    temp1 = ismember(ABDIr(i).Inputs,ALX.cellIDs');
    if sum(temp1)~=0 ;
        ALX.ABDIrpathLength(i,:) = histcounts(ABDIr(i).PathLength(temp1)'./(max(Pvec_tree(ABDIr(i).Tree{1}))),0:0.1:1,'Normalization','probability');
    else
        ALX.ABDIrpathLength(i,:) = repmat(NaN,1,10);
    end
    
end


for i = 1:numel(ABDIc)
    temp1 = ismember(ABDIc(i).Inputs,ALX.cellIDs');
    if sum(temp1)~=0 ;
        ALX.ABDIcpathLength(i,:) = histcounts(ABDIc(i).PathLength(temp1)'./(max(Pvec_tree(ABDIc(i).Tree{1}))),0:0.1:1,'Normalization','probability');
    else
        ALX.ABDIcpathLength(i,:) = repmat(NaN,1,10);
    end
end


subplot(4,4,1)

ALX.ABDpathLength = vertcat(ALX.ABDrpathLength,ALX.ABDcpathLength);
ALX.ABDipathLength = vertcat(ALX.ABDIrpathLength,ALX.ABDIcpathLength);

ALX.ABDpathLength = ALX.ABDpathLength(ALX.ABDpathLength~=0);
ALX.ABDipathLength = ALX.ABDipathLength(ALX.ABDipathLength~=0);

errorbar(nanmean(vertcat(ALX.ABDrpathLength,ALX.ABDcpathLength)),nanstd([ALX.ABDrpathLength;ALX.ABDcpathLength])./sqrt(29),...
    '-o','color',ALXcolor,'LineWidth',2,'MarkerFaceColor','w');
hold on;
errorbar(nanmean(vertcat(ALX.ABDIrpathLength,ALX.ABDIcpathLength)),nanstd([ALX.ABDIrpathLength;ALX.ABDIcpathLength])./sqrt(29),...
    '-o','color',ALXcolor,'LineWidth',2,'MarkerFaceColor',ALXcolor);
axis square;
legend({'onto ABD (M)', 'onto ABDi (I)'},'Location','bestoutside');
set(gca,'XTickLabels',[0,0.5,1]);
xlabel('Norm. pathlength');
ylabel('Norm. count');
box off;
offsetAxes(gca);




%%

allContraAxons.cellID = vertcat(Int.Contra);
allContraAxons.cellID = unique(allContraAxons.cellID);
allContraAxons.rhombomeres = isRhombomere(allContraAxons.cellID);
allContraAxons.r6r7 = [allContraAxons.cellID(logical(allContraAxons.rhombomeres.r6));...
    allContraAxons.cellID(logical(allContraAxons.rhombomeres.r7))];

allContraAxons.r6r7MotorDist = isMotor(allContraAxons.r6r7,df);
allContraAxons.r6r7isMotor = allContraAxons.cellID(sum(allContraAxons.r6r7MotorDist(:,2:end),2)>0);

allContraAxons.r6r7NormCount = (sum(allContraAxons.r6r7MotorDist(:,2:3),2) - sum(allContraAxons.r6r7MotorDist(:,4:5),2)) ./ ...
    (sum(allContraAxons.r6r7MotorDist(:,2:3),2) + sum(allContraAxons.r6r7MotorDist(:,4:5),2));

figure;
transform_swc_AV(allContraAxons.r6r7isMotor,'k',[],true,false);

%
allDBXcontra.cellID = vertcat(Int(45:75).Contra);
allDBXcontra.cellID = unique(allDBXcontra.cellID);
allDBXcontra.cellID =  allDBXcontra.cellID(isExistReRoot(allDBXcontra.cellID));
allDBXcontra.isMotor = isMotor(allDBXcontra.cellID ,df);
allDBXcontra.cellID = allDBXcontra.cellID(sum(allDBXcontra.isMotor(:,2:end),2)>0);

% manual list

allDBXcontra.CellID = [77865,76776,77773,78651,81356,77444,78631,77370,77811,81682,76843,77852,76774,77819,...
80246,77768,80248,81397,81385,76782,81840,78653,81681,76666,81597,78677,78227,78671,...
76693,80219,77358,81106,78687,77830,77771,80575,81138,80203,80212,77774,78655];

vestibularCells = unique(vertcat(Int.Vestibular));

for i = 1:numel(vestibularCells)
    Vest(i) = InputsByClass(vestibularCells(i),df,1);
end

[allDBXcontra.VestibularPathLength,allDBXcontra.VestibularGradient] = getABDgradient(Vest,allDBXcontra.CellID',true,false);


%% ALX Pairwise distances
for i = 1:numel(ALX.cellIDs)
    for j = 1:numel(ALX.cellIDs)
        
        indexi = find(ismember([Int.cellID]',ALX.cellIDs(i)));
        indexj = find(ismember([Int.cellID]',ALX.cellIDs(j)));
        
        ALX.pariwiseDist(i,j) = sqrt(sum((ALX.Origin(i,:) - ALX.Origin(j,:)).^2));
        commonSaccadic = intersect(Int(indexi).Saccadic,Int(indexj).Saccadic);
        commonVestibular = intersect(Int(indexi).Vestibular,Int(indexj).Vestibular);
        commonIntgrators = intersect(Int(indexi).Integrator,Int(indexj).Integrator);
        commonContra = intersect(Int(indexi).Contra,Int(indexj).Contra);
        
        ALX.commonSaccadic(i,j) = length(commonSaccadic);
        ALX.commonVestibular(i,j) = length(commonVestibular);
        ALX.commonIntegrators(i,j) = length(commonIntgrators);
        ALX.commonContra(i,j) = length(commonContra);
        
        clear commonSaccadic;
        clear commonVestibular;
        clear commonIntgrators;
        clear commonContra;
    end
end

ALX.commonSaccadic(ALX.commonSaccadic ==0) = NaN;
ALX.commonVestibular(ALX.commonVestibular ==0) = NaN;
ALX.commonIntegrators(ALX.commonIntegrators ==0) = NaN;
ALX.commonContra(ALX.commonContra ==0) = NaN;
ALX.pariwiseDist(ALX.pariwiseDist ==0) = NaN;
[h,e,ind] = histcounts(ALX.pariwiseDist(:),0:10:120);

figure;
subplot(4,4,1)
for i = 1:max(ind)
    ALXSaccadic(i,:) = [nanmean(ALX.pariwiseDist(i==ind)),nanmean(ALX.commonSaccadic(i==ind))];
    errorbar(nanmean(ALX.pariwiseDist(i==ind)),nanmean(ALX.commonSaccadic(i==ind)),...
        nanstd(ALX.commonSaccadic(i==ind))./sqrt(numel(ALX.commonSaccadic(i==ind))),'vertical',...
        'o','MarkerFaceColor','w','MarkerEdgeColor',colors(1,:),'MarkerSize',4,'Color',colors(1,:),'LineWidth',2);
    hold on
end


axis square
xlabel('Pairwise distance (\mum)');
ylabel('Common axons');
title('Saccadic')
box off;

subplot(4,4,2)
for i = 1:max(ind)
    ALXVestibular(i,:) = [nanmean(ALX.pariwiseDist(i==ind)),nanmean(ALX.commonVestibular(i==ind))];
    errorbar(nanmean(ALX.pariwiseDist(i==ind)),nanmean(ALX.commonVestibular(i==ind)),...
        nanstd(ALX.commonVestibular(i==ind))./sqrt(numel(ALX.commonVestibular(i==ind))),'vertical',...
        'o','MarkerFaceColor','w','MarkerEdgeColor',colors(1,:),'MarkerSize',4,'Color',colors(1,:),'LineWidth',2);
    hold on
    
end

axis square
xlabel('Pairwise distance (\mum)');
ylabel('Common axons');
title('Vestibular')
box off;

subplot(4,4,3)
for i = 1:max(ind)
    ALXIntegrator(i,:) = [nanmean(ALX.pariwiseDist(i==ind)),nanmean(ALX.commonIntegrators(i==ind))];
    errorbar(nanmean(ALX.pariwiseDist(i==ind)),nanmean(ALX.commonIntegrators(i==ind)),...
        nanstd(ALX.commonIntegrators(i==ind))./sqrt(numel(ALX.commonIntegrators(i==ind))),'vertical',...
        'o','MarkerFaceColor','w','MarkerEdgeColor',colors(1,:),'MarkerSize',4,'Color',colors(1,:),'LineWidth',2);
    hold on
    
end

axis square
xlabel('Pairwise distance (\mum)');
ylabel('Common axons');
title('Integrators')
box off;


subplot(4,4,4)
for i = 1:max(ind)
    ALXContra(i,:) = [nanmean(ALX.pariwiseDist(i==ind)),nanmean(ALX.commonContra(i==ind))];
    errorbar(nanmean(ALX.pariwiseDist(i==ind)),nanmean(ALX.commonContra(i==ind)),...
        nanstd(ALX.commonContra(i==ind))./sqrt(numel(ALX.commonContra(i==ind))),'vertical',...
        'o','MarkerFaceColor','w','MarkerEdgeColor',colors(1,:),'MarkerSize',4,'Color',colors(1,:),'LineWidth',2);
    hold on
    
end

axis square
xlabel('Pairwise distance (\mum)');
ylabel('Common axons');
title('Contra')
box off;

%% DBX Pairwise distances
for i = 1:numel(DBX.cellIDs)
    for j = 1:numel(DBX.cellIDs)
        
        indexi = find(ismember([Int.cellID]',DBX.cellIDs(i)));
        indexj = find(ismember([Int.cellID]',DBX.cellIDs(j)));
        
        DBX.pariwiseDist(i,j) = sqrt(sum((DBX.Origin(i,:) - DBX.Origin(j,:)).^2));
        commonSaccadic = intersect(Int(indexi).Saccadic,Int(indexj).Saccadic);
        commonVestibular = intersect(Int(indexi).Vestibular,Int(indexj).Vestibular);
        commonIntgrators = intersect(Int(indexi).Integrator,Int(indexj).Integrator);
        commonContra = intersect(Int(indexi).Contra,Int(indexj).Contra);
        
        DBX.commonSaccadic(i,j) = length(commonSaccadic);
        DBX.commonVestibular(i,j) = length(commonVestibular);
        DBX.commonIntegrators(i,j) = length(commonIntgrators);
        DBX.commonContra(i,j) = length(commonContra);
        
        clear commonSaccadic;
        clear commonVestibular;
        clear commonIntgrators;
        clear commonContra;
    end
end

DBX.commonSaccadic(DBX.commonSaccadic ==0) = NaN;
DBX.commonVestibular(DBX.commonVestibular ==0) = NaN;
DBX.commonIntegrators(DBX.commonIntegrators ==0) = NaN;
DBX.commonContra(DBX.commonContra ==0) = NaN;
DBX.pariwiseDist(DBX.pariwiseDist ==0) = NaN;
[h,e,ind] = histcounts(DBX.pariwiseDist(:),0:10:120);

hold on;
subplot(4,4,1)
for i = 1:max(ind)
    DBXSaccadic(i,:) = [nanmean(DBX.pariwiseDist(i==ind)),nanmean(DBX.commonSaccadic(i==ind))];
    errorbar(nanmean(DBX.pariwiseDist(i==ind)),nanmean(DBX.commonSaccadic(i==ind)),...
        nanstd(DBX.commonSaccadic(i==ind))./sqrt(numel(DBX.commonSaccadic(i==ind))),'vertical',...
        'o','MarkerFaceColor','w','MarkerEdgeColor',colors(2,:),'MarkerSize',4,'Color',colors(2,:),'LineWidth',2);
    hold on
end


axis square
xlabel('Pairwise distance (\mum)');
ylabel('Common axons');
title('Saccadic')
box off;

subplot(4,4,2)
for i = 1:max(ind)
    DBXVestibular(i,:) = [nanmean(DBX.pariwiseDist(i==ind)),nanmean(DBX.commonVestibular(i==ind))];
    errorbar(nanmean(DBX.pariwiseDist(i==ind)),nanmean(DBX.commonVestibular(i==ind)),...
        nanstd(DBX.commonVestibular(i==ind))./sqrt(numel(DBX.commonVestibular(i==ind))),'vertical',...
        'o','MarkerFaceColor','w','MarkerEdgeColor',colors(2,:),'MarkerSize',4,'Color',colors(2,:),'LineWidth',2);
    hold on
    
end

axis square
xlabel('Pairwise distance (\mum)');
ylabel('Common axons');
title('Vestibular')
box off;

subplot(4,4,3)
for i = 1:max(ind)
    DBXIntegrator(i,:) = [nanmean(DBX.pariwiseDist(i==ind)),nanmean(DBX.commonIntegrators(i==ind))];
    errorbar(nanmean(DBX.pariwiseDist(i==ind)),nanmean(DBX.commonIntegrators(i==ind)),...
        nanstd(DBX.commonIntegrators(i==ind))./sqrt(numel(DBX.commonIntegrators(i==ind))),'vertical',...
        'o','MarkerFaceColor','w','MarkerEdgeColor',colors(2,:),'MarkerSize',4,'Color',colors(2,:),'LineWidth',2);
    hold on
    
end

axis square
xlabel('Pairwise distance (\mum)');
ylabel('Common axons');
title('Integrators')
box off;


subplot(4,4,4)
for i = 1:max(ind)
    DBXContra(i,:) = [nanmean(DBX.pariwiseDist(i==ind)),nanmean(DBX.commonContra(i==ind))];
    errorbar(nanmean(DBX.pariwiseDist(i==ind)),nanmean(DBX.commonContra(i==ind)),...
        nanstd(DBX.commonContra(i==ind))./sqrt(numel(DBX.commonContra(i==ind))),'vertical',...
        'o','MarkerFaceColor','w','MarkerEdgeColor',colors(2,:),'MarkerSize',4,'Color',colors(2,:),'LineWidth',2);
    hold on
    
end

axis square
xlabel('Pairwise distance (\mum)');
ylabel('Common axons');
title('Contra')
box off;

%%

subplot(4,4,1)

showfit(ezfit(ALXSaccadic(:,1),ALXSaccadic(:,2),'power'),'fitcolor',colors(1,:),'dispeqboxmode','off');
showfit(ezfit(DBXSaccadic(:,1),DBXSaccadic(:,2),'power'),'fitcolor',colors(2,:),'dispeqboxmode','off');

subplot(4,4,2)
showfit(ezfit(ALXVestibular(:,1),ALXVestibular(:,2),'power'),'fitcolor',colors(1,:),'dispeqboxmode','off');
showfit(ezfit(DBXVestibular(:,1),DBXVestibular(:,2),'power'),'fitcolor',colors(2,:),'dispeqboxmode','off');


subplot(4,4,3)
showfit(ezfit(ALXIntegrator(:,1),ALXIntegrator(:,2),'power'),'fitcolor',colors(1,:),'dispeqboxmode','off');
showfit(ezfit(DBXIntegrator(1:end-2,1),DBXIntegrator(1:end-2,2),'power'),'fitcolor',colors(2,:),'dispeqboxmode','off');

subplot(4,4,4)
showfit(ezfit(ALXContra(:,1),ALXContra(:,2),'power'),'fitcolor',colors(1,:),'dispeqboxmode','off');
showfit(ezfit(DBXContra(1:end-2,1),DBXContra(1:end-2,2),'power'),'fitcolor',colors(2,:),'dispeqboxmode','off');


%% Integrator funcitons

load STA.mat
load FunctionsCellIDs.mat
load('MelanieDBXCells.mat');
load('MelanieALXCells.mat');
load('IntegratorsSaccABD.mat');
load('IntegratorsSaccABDi.mat');

t = [-2:0.05:7];


Firing = [STAall, DBX_vglut_neg/100 , DBX_vglut/100, ALX_neg/100, ALX_pos/100 ]; % divide by 100 to convert back from %

% consider only the first 5 components of the trace
[A,B,C,D,E,F] = pca(Firing);
% sprintf('first %d components capture %f of the data',2, sum(E(1:2)))
STApc = B(:,1:2)*A(:,1:2)'+ F;
STApc = STApc(:,1:22);

for i = 1:size(STApc,2)
    STAdeconv(:,i) = DeconvCa(STApc(:,i));
end

FunctionalIntegratorsSaccABD = intersect(IntegratorsSaccABD,functionalCellIDs_new);
FunctionalIntegratorsSaccABDi = intersect(IntegratorsSaccABDi,functionalCellIDs_new);

confirmedIntABD = find(ismember(functionalCellIDs_new,FunctionalIntegratorsSaccABD));
confirmedIntABDi = find(ismember(functionalCellIDs_new,FunctionalIntegratorsSaccABDi));

figure;
subplot(4,4,1)
plot(t,STApc(:,confirmedIntABD)./max(STApc(:,confirmedIntABD)),'color',leadColor);
box off;
axis square;

subplot(4,4,3)
plot(t,STApc(:,confirmedIntABDi)./max(STApc(:,confirmedIntABDi)),'color',lagColor);
box off;
axis square;

subplot(4,4,5)
shadedErrorBar(t,mean(STApc(:,confirmedIntABD),2),std(STApc(:,confirmedIntABD),[],2),'lineprops',{'Color',leadColor});
shadedErrorBar(t,mean(STApc(:,confirmedIntABDi),2),std(STApc(:,confirmedIntABDi),[],2),'lineprops',{'Color',lagColor});
axis square;
box off;

subplot(4,4,7)
shadedErrorBar(t,mean(STAdeconv(:,confirmedIntABD),2),std(STAdeconv(:,confirmedIntABD),[],2),'lineprops',{'Color',leadColor});
shadedErrorBar(t,mean(STAdeconv(:,confirmedIntABDi),2),std(STAdeconv(:,confirmedIntABDi),[],2),'lineprops',{'Color',lagColor});
axis square;
box off;


%%

[~,FunctionalIntegratorsSaccABDorder] = intersect(IntegratorNeurons,FunctionalIntegratorsSaccABD);
[~,FunctionalIntegratorsSaccABDiorder] = intersect(IntegratorNeurons,FunctionalIntegratorsSaccABDi);

figure;
for i = 1:numel(FunctionalIntegratorsSaccABDorder)
    subplot(4,4,i)
    plot(t,STApc(:,confirmedIntABD(i))./max(STApc(:,confirmedIntABD(i))),'color',leadColor);
    int = sprintf('integrators:%d \n Saccadics:%d',numel(Int(FunctionalIntegratorsSaccABDorder(i)).Integrator),...
        numel(Int(FunctionalIntegratorsSaccABDorder(i)).Saccadic));
    text(6.5,1,int);
end

for j = 1:numel(FunctionalIntegratorsSaccABDiorder)
    subplot(4,4,i+j)
    plot(t,STApc(:,confirmedIntABDi(j))./max(STApc(:,confirmedIntABDi(j))),'color',lagColor);
    int = sprintf('integrators:%d \n Saccadics:%d',numel(Int(FunctionalIntegratorsSaccABDiorder(j)).Integrator), ...
        numel(Int(FunctionalIntegratorsSaccABDiorder(i)).Saccadic));
    text(6.5,1,int);
end

figure;

AllFunctionalIntegrators = [FunctionalIntegratorsSaccABD;FunctionalIntegratorsSaccABDi];
AllFunctionalIntegratorsOrder = [FunctionalIntegratorsSaccABDorder;FunctionalIntegratorsSaccABDiorder];
AllConfirmedInt = [confirmedIntABD';confirmedIntABDi'];
newT = 0:0.05:0.05*(181-41);

for i = 1:numel(AllFunctionalIntegratorsOrder)
    NumberOfIntegratorSynapses(i) = numel(Int(AllFunctionalIntegratorsOrder(i)).Integrator);
    F = ezfit(newT, STApc(41:end,i),'a*exp(-x/t)')
    Tau(i) = F.m(2);
end
Tau(Tau>100)=100;
Tau(Tau<1)=1;
[~,NumberOfIntegratorSynapsesSorted] = sort(NumberOfIntegratorSynapses);

for i = 1:numel(NumberOfIntegratorSynapsesSorted)
    plot(t,0.5*i+ STApc(:,AllConfirmedInt(NumberOfIntegratorSynapsesSorted(i)))./max(STApc(:,AllConfirmedInt(NumberOfIntegratorSynapsesSorted(i)))));
    hold on;
end



