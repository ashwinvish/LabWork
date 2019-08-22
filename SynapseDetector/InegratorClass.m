% Classes by anatomy
%clear;
addpath(genpath('/Users/ashwin/Documents/'));

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

ALX.cellIDs = [confirmedALX,putativeALX];
ALX.cellIDs = ALX.cellIDs(isExistReRoot(ALX.cellIDs));

confirmedDBX = [76199 76200 76182 76183 76185 76186 76189 76191 76188 ];
putativeDBX = [77582 77588 77595 77597 77598 77599 77605 79040 76399 76828 ...
    76829 78435 76826 76874 76289 76542 76832 76836 76838 76877 76883 76884 ...
    78400 78429 76383 76867 77594 76830 76876 78440 78438 78447 79045 78434 ...
    78450 76878 77683 77450 78448 ];

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
     Int(i) = InputsByClass(IntegratorNeurons(i),df);
     % locate and remove self touches among the integrators
     ind = find(ismember(Int(i).Integrator,IntegratorNeurons(i))); 
     Int(i).Integrator(ind) = [];
     clear ind;
     end  
end
 
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

for i = 1:numel(DBX.cellIDs)
    indexi = find(ismember([Int.cellID]',DBX.cellIDs(i)));
    if ~isempty(indexi)
    DBX.SaccadicNumbers(i) = length(Int(indexi).Saccadic);
    DBX.VestibularNumbers(i) = length(Int(indexi).Vestibular);
    DBX.IntegratorNumbers(i) = length(Int(indexi).Integrator);
    DBX.ContraNumbers(i) = length(Int(indexi).Contra);
    DBX.RestNumbers(i) = length(Int(indexi).EverythingElse);
    DBX.Origin(i,:) = Int(indexi).Origin
    end
end

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
xlabel('(A-Ai)/(A+Ai)');
axis square;
box off;

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

%%

allContraAxons.cellID = vertcat(Int.Contra);
allContraAxons.cellID = unique(allContraAxons.cellID);
allContraAxons.rhombomeres = isRhombomere(allContraAxons.cellID);
allContraAxons.r6r7 = [allContraAxons.cellID(logical(allContraAxons.rhombomeres.r6));...
                        allContraAxons.cellID(logical(allContraAxons.rhombomeres.r7))];

allContraAxons.r6r7MotorDist = isMotor(allContraAxons.r6r7,df);

allContraAxons.r6r7NormCount = (sum(allContraAxons.r6r7MotorDist(:,2:3),2) - sum(allContraAxons.r6r7MotorDist(:,4:5),2)) ./ ...
                            (sum(allContraAxons.r6r7MotorDist(:,2:3),2) + sum(allContraAxons.r6r7MotorDist(:,4:5),2));



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



