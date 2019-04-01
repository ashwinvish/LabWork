% Burst vs Tonic Connectomes

clear;
addpath(genpath('/Users/ashwin/Documents/'));

colors = cbrewer('qual','Dark2',10);
startup

if ismac
    addpath(genpath('/Users/ashwin/Documents/'));
    df = readtable('/Users/ashwin/Documents/SynapseDetector/11252018.csv');
else
    addpath(genpath('/usr/people/ashwinv/seungmount/research/Ashwin/Scripts'));
    df = readtable('/usr/people/ashwinv/seungmount/research/Ashwin/SynapseDetector/11252018.csv');
end

%load('FiringRates.mat');
load('TAU.mat');
load('STA.mat');
load('AllCells.mat');
load('MelanieDBXCells.mat');
load('MelanieALXCells.mat');



DbxCells = [76182,76183,76185,76186,76188,76189,76191,76199,76200];
OriginalCellOrderDBX = [8,9,10,11,21,12,15,2,3];
DbxTimeConstants = TAU(OriginalCellOrderDBX);

t = [-2:0.05:7]; % time period in sec

Firing = [STAall(:,OriginalCellOrderDBX) , DBX_vglut_neg/100 , DBX_vglut/100]; % divide by 100 to convert back from %
DBXpop = [DBX_vglut_neg/100 , DBX_vglut/100]; % melanies data is reported as %
CellOrder = [ones(1,size(DbxCells,2)), 2* ones(1,size(DBX_vglut_neg,2)), 3* ones(1,size(DBX_vglut,2))];

% consider only the first 5 components of the trace
[A,B,C,D,E,F] = pca(Firing);
sprintf('first %d components capture %f of the data',5, sum(E(1:5)))
Firing = B(:,1:5)*A(:,1:5)'+ F;
normFiring = Firing(:,1:9);


temp1 = colorcet('CBTL1','N',5); %  linear-tritanopic_krjcw_5-98_c46_n256
temp2 = colorcet('CBL2','N',5);  %  linear-protanopic-deuteranopic_kbw_5-98_c40_n256

leadColor = temp1(3,:); % dark red, CB safe
lagColor = temp2(3,:); % dark blue, CB safe

clear temp1;
clear temp2;

PartnerColors = colorcet('CBTD1','N',10);
lightRed = PartnerColors(10,:);
lightBlue = PartnerColors(1,:);


molOrder = [find(CellOrder ==1)'; find(CellOrder ==2)'; find(CellOrder ==3)'];
eva = evalclusters(Firing','kmeans','silhouette','Distance','correlation','KList',1:6);

order = [find(eva.OptimalY==1);find(eva.OptimalY ==2)];

if max(mean(Firing(1:50,eva.OptimalY==1),2)) > max(mean(Firing(1:50,eva.OptimalY==2),2))
    leadCorrIDs = find(eva.OptimalY==1);
    lagCorrIDs = find(eva.OptimalY==2);
else
    leadCorrIDs = find(eva.OptimalY==2);
    lagCorrIDs = find(eva.OptimalY==1);
end


leadNeurons = DbxCells(leadCorrIDs(leadCorrIDs<10));
lagNeurons  = DbxCells(lagCorrIDs(lagCorrIDs<10));


% lead Neurons
for i = 1:numel(leadNeurons)
    Lead(i) = InputsByClass(leadNeurons(i),df);
end

% lag neurons
for i = 1:numel(lagNeurons)
    Lag(i)  = InputsByClass(lagNeurons(i),df);
end

%% Burst/Tonic Connectome


load ConnMatrixPre.mat;
load AllCells.mat;


SaccadicNeurons = vertcat(Lead.Saccadic,Lag.Saccadic);
SaccadicNeurons = unique(SaccadicNeurons(SaccadicNeurons<1e5),'stable');
VestibularNeuron = vertcat(Lead.Vestibular,Lag.Vestibular);
VestibularNeuron = unique(VestibularNeuron(VestibularNeuron<1e5),'stable');
IntegratorNeurons = vertcat(Lead.Integrator, Lag.Integrator);
IntegratorNeurons = unique(IntegratorNeurons(IntegratorNeurons<1e5),'stable');
RemainingNeurons = vertcat(Lead.EverythingElse, Lag.EverythingElse);
RemainingNeurons = unique(RemainingNeurons(RemainingNeurons<1e5),'stable');
ContraNeurons = vertcat(Lead.Contra, Lag.Contra);
ContraNeurons = unique(ContraNeurons(ContraNeurons<1e5),'stable');

AllInputs = [SaccadicNeurons;VestibularNeuron;IntegratorNeurons;RemainingNeurons;ContraNeurons];
AllInputs = unique(AllInputs(AllInputs<1e5),'stable');

LeadInputs = vertcat(Lead.Inputs);
LeadInputs = unique(LeadInputs(LeadInputs<1e5),'stable');

LagInputs = vertcat(Lag.Inputs);
LagInputs = unique(LagInputs(LagInputs<1e5),'stable');

[~,~,allInputLocation] = intersect(AllInputs,AllCells,'stable');

allInputConnectome = zeros(size(allInputLocation,1));

for i = 1:size(allInputLocation,1)
    allInputConnectome(i,:) = ConnMatrixPre(allInputLocation(i),allInputLocation);
end

[~,~,leadA] = intersect(LeadInputs,AllInputs,'stable');

BurstConnetome = zeros(size(AllInputs,1));

for i = 1:size(leadA,1)
    BurstConnetome(leadA(i),leadA) = allInputConnectome(leadA(i),leadA);
end


[~,~,lagA] = intersect(LagInputs,AllInputs,'stable');


TonicConnectome = zeros(size(AllInputs,1));

for i = 1:size(lagA,1)
    TonicConnectome(lagA(i),lagA) = allInputConnectome(lagA(i),lagA);
end

%%
SaccadicBloc = [1,find(AllInputs == SaccadicNeurons(end))];
SaccadicBloc = SaccadicBloc+0.5;
VestibularBloc = [find(AllInputs == SaccadicNeurons(end))+1 ,find(AllInputs == VestibularNeuron(end))];
VestibularBloc = VestibularBloc +0.5;
IntegratorBloc = [find(AllInputs == VestibularNeuron(end))+1 , find(AllInputs == IntegratorNeurons(end))];
IntegratorBloc = IntegratorBloc+0.5;
EverythingElseBloc = [find(AllInputs == IntegratorNeurons(end))+1 , find(AllInputs == RemainingNeurons(end))];
EverythingElseBloc = EverythingElseBloc+0.5;
ContraBloc = [find(AllInputs == RemainingNeurons(end))+1 , find(AllInputs == ContraNeurons(end))];

subplot(1,2,1);
maxRed = max(BurstConnetome(:));
%redMap = cbrewer('seq','Reds',maxRed);
redMap = colorcet('L4','N',maxRed);
h1 = cspy(BurstConnetome,'ColorMap',redMap,'Level',maxRed,'MarkerSize',20);
%imagesc(BurstConnetome);
colormap(redMap);

line([SaccadicBloc(2),SaccadicBloc(2)],[0,size(BurstConnetome,1)],'Color','k');
line([0,size(BurstConnetome,1)],[SaccadicBloc(2),SaccadicBloc(2)],'Color','k');

line([VestibularBloc(2),VestibularBloc(2)],[0,size(BurstConnetome,1)],'Color','k');
line([0,size(BurstConnetome,1)],[VestibularBloc(2),VestibularBloc(2)],'Color','k');

line([IntegratorBloc(2),IntegratorBloc(2)],[0,size(BurstConnetome,1)],'Color','k');
line([0,size(BurstConnetome,1)],[IntegratorBloc(2),IntegratorBloc(2)],'Color','k');


line([EverythingElseBloc(2),EverythingElseBloc(2)],[0,size(BurstConnetome,1)],'Color','k');
line([0,size(BurstConnetome,1)],[EverythingElseBloc(2),EverythingElseBloc(2)],'Color','k');

title('Burst type');
colorbar(h1);
box on;
set(h1,'XTick',[],'YTick',[]);

subplot(1,2,2);
maxBlue = max(TonicConnectome(:));
blueMap = colorcet('L6','N',maxBlue);
h2 = cspy(TonicConnectome,'ColorMap',blueMap,'Level',maxBlue,'MarkerSize',20);
%imagesc(TonicConnectome);
colormap(h2,blueMap);

line([SaccadicBloc(2),SaccadicBloc(2)],[0,size(BurstConnetome,1)],'Color','k');
line([0,size(BurstConnetome,1)],[SaccadicBloc(2),SaccadicBloc(2)],'Color','k');

line([VestibularBloc(2),VestibularBloc(2)],[0,size(BurstConnetome,1)],'Color','k');
line([0,size(BurstConnetome,1)],[VestibularBloc(2),VestibularBloc(2)],'Color','k');

line([IntegratorBloc(2),IntegratorBloc(2)],[0,size(BurstConnetome,1)],'Color','k');
line([0,size(BurstConnetome,1)],[IntegratorBloc(2),IntegratorBloc(2)],'Color','k');


line([EverythingElseBloc(2),EverythingElseBloc(2)],[0,size(BurstConnetome,1)],'Color','k');
line([0,size(BurstConnetome,1)],[EverythingElseBloc(2),EverythingElseBloc(2)],'Color','k');

title('Tonic type');
colorbar(h2);
box on;
set(h2,'XTick',[],'YTick',[]);
colormap(h1,redMap);
