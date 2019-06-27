
clear MatOrder
clear MatIndex
clear connMat

load ABDr.mat
load ABDc.mat
load ABDIr.mat
load ABDIc.mat

load AllCells.mat
load ConnMatrixPre.mat

load ABDVestPop.mat
load ABDiVestPop.mat
load ABDContraPop.mat
load ABDiContrapop.mat

addpath(genpath('/Users/ashwin/Documents/'));

colors = cbrewer('qual','Dark2',10);

temp1 = colorcet('CBTL1','N',5); %  linear-tritanopic_krjcw_5-98_c46_n256
temp2 = colorcet('CBL2','N',5);  %  linear-protanopic-deuteranopic_kbw_5-98_c40_n256

leadColor = temp1(3,:); % dark red, CB safe
lagColor = temp2(3,:); % dark blue, CB safe

PartnerColors = colorcet('CBTD1','N',10);
lightRed = PartnerColors(10,:);
lightBlue = PartnerColors(1,:);


if ismac
    addpath(genpath('/Users/ashwin/Documents/'));
    df = readtable('/Users/ashwin/Documents/SynapseDetector/04152019.csv');
else
    addpath(genpath('/usr/people/ashwinv/seungmount/research/Ashwin/Scripts'));
    df = readtable('/usr/people/ashwinv/seungmount/research/Ashwin/SynapseDetector/11252018.csv');
end

%%

ABDr_uniqueInputs = unique(vertcat(ABDr.Inputs));
ABDr_uniqueInputs = ABDr_uniqueInputs(ABDr_uniqueInputs<1e5);


ABDc_uniqueInputs = unique(vertcat(ABDc.Inputs));
ABDc_uniqueInputs = ABDc_uniqueInputs(ABDc_uniqueInputs<1e5);


ABDIr_uniqueInputs = unique(vertcat(ABDIr.Inputs));
ABDIr_uniqueInputs = ABDIr_uniqueInputs(ABDIr_uniqueInputs<1e5);


ABDIc_uniqueInputs = unique(vertcat(ABDIc.Inputs));
ABDIc_uniqueInputs = ABDIc_uniqueInputs(ABDIc_uniqueInputs<1e5);


All.Saccadic = [unique(vertcat(ABDr.Saccadic));unique(vertcat(ABDc.Saccadic));...
    unique(vertcat(ABDIr.Saccadic));unique(vertcat(ABDIc.Saccadic))];
All.Vestibular = [unique(vertcat(ABDr.Vestibular));unique(vertcat(ABDc.Vestibular));...
    unique(vertcat(ABDIr.Vestibular));unique(vertcat(ABDIc.Vestibular))];
All.Integrator = [unique(vertcat(ABDr.Integrator));unique(vertcat(ABDc.Integrator));...
    unique(vertcat(ABDIr.Integrator));unique(vertcat(ABDIc.Integrator))];
All.Contra = [unique(vertcat(ABDr.Contra));unique(vertcat(ABDc.Contra));...
    unique(vertcat(ABDIr.Contra));unique(vertcat(ABDIc.Contra))];
All.Rest = [unique(vertcat(ABDr.EverythingElse));unique(vertcat(ABDc.EverythingElse));...
    unique(vertcat(ABDIr.EverythingElse));unique(vertcat(ABDIc.EverythingElse))];
All.ABD = [vertcat(ABDr.cellID);vertcat(ABDc.cellID);vertcat(ABDIr.cellID);vertcat(ABDIc.cellID)]

All.Inputs = unique([ABDr_uniqueInputs ;ABDc_uniqueInputs; ABDIr_uniqueInputs ;ABDIc_uniqueInputs]);


All.SaccadicIndex = find(ismember(All.Inputs,All.Saccadic));
All.VestibularIndex = find(ismember(All.Inputs,All.Vestibular));
All.IntegratorIndex = find(ismember(All.Inputs,All.Integrator));
All.ContraIndex = find(ismember(All.Inputs,All.Contra));
All.RestIndex = find(ismember(All.Inputs,All.Rest));


All.MotorCounts = isMotor(All.Inputs,df);
All.ABDcounts = sum(All.MotorCounts(:,2:3),2);
All.ABDicounts = sum(All.MotorCounts(:,4:5),2);
All.normalizedMotorCounts = (All.ABDcounts-All.ABDicounts)./ (All.ABDcounts+All.ABDicounts);
[~,All.normalizedMotorCountsSortedIndex] = sort(All.normalizedMotorCounts);


figure;
subplot(4,4,1)
histogram(All.normalizedMotorCounts,-1:0.1:1,'FaceColor','k');
axis square;
box off;
title('All');

subplot(4,4,3)
histogram(All.normalizedMotorCounts(All.SaccadicIndex),-1:0.1:1,'FaceColor','k');
axis square;
box off;
title('Saccadic');

subplot(4,4,5)
histogram(All.normalizedMotorCounts(All.VestibularIndex),-1:0.1:1,'FaceColor','k');
axis square;
box off;
title('Vestibular');


subplot(4,4,7)
histogram(All.normalizedMotorCounts(All.IntegratorIndex),-1:0.1:1,'FaceColor','k');
axis square;
box off;
title('Integrator');


subplot(4,4,9)
histogram(All.normalizedMotorCounts(All.ContraIndex),-1:0.1:1,'FaceColor','k');
axis square;
box off;
title('Contra');


subplot(4,4,11)
histogram(All.normalizedMotorCounts(All.RestIndex),-1:0.1:1,'FaceColor','k');
axis square;
box off;
title('Rest');

%%

All.SaccadicABD = All.Inputs(intersect(find(All.normalizedMotorCounts>0),All.SaccadicIndex));
All.SaccadicABDi = All.Inputs(intersect(find(All.normalizedMotorCounts<0),All.SaccadicIndex));

All.VestibularABD = All.Inputs(intersect(find(All.normalizedMotorCounts>0),All.VestibularIndex));
All.VestibularABDi = All.Inputs(intersect(find(All.normalizedMotorCounts<0),All.VestibularIndex));

All.IntegratorABD = All.Inputs(intersect(find(All.normalizedMotorCounts>0),All.IntegratorIndex));
All.IntegratorABDi = All.Inputs(intersect(find(All.normalizedMotorCounts<0),All.IntegratorIndex));

All.ContraABD = All.Inputs(intersect(find(All.normalizedMotorCounts>0),All.ContraIndex));
All.ContraABDi = All.Inputs(intersect(find(All.normalizedMotorCounts<0),All.ContraIndex));

All.RestABD = All.Inputs(intersect(find(All.normalizedMotorCounts>0),All.RestIndex));
All.RestABDi = All.Inputs(intersect(find(All.normalizedMotorCounts<0),All.RestIndex));

% remove ABD contamination from Rest
All.RestABD = All.RestABD(~ismember(All.RestABD,All.ABD));
All.RestABDi = All.RestABDi(~ismember(All.RestABDi,All.ABD));

ABDRestPop = All.RestABD;
ABDiRestPop = All.RestABDi;

save('ABDRestPop.mat','ABDRestPop');
save('ABDiRestPop.mat','ABDiRestPop');
%%

ALL.ABDonly = All.Inputs(All.normalizedMotorCounts==1);
ALL.ABDionly = All.Inputs(All.normalizedMotorCounts==-1);

subplot(1,2,1)
transform_swc_AV(ALL.ABDonly(isExistReRoot(ALL.ABDonly)),leadColor,[],true,false);
subplot(1,2,2)
transform_swc_AV(ALL.ABDionly(isExistReRoot(ALL.ABDionly)),lagColor,[],true,false);

%% determine rhombomere organization

All.ContraABDRhombomere = isRhombomere(All.ContraABD);
All.ContraABDiRhombomere = isRhombomere(All.ContraABDi);

All.ContraABDOrdered = [All.ContraABDRhombomere.cellID(All.ContraABDRhombomere.r3==1);...
                        All.ContraABDRhombomere.cellID(All.ContraABDRhombomere.r4==1);...
                        All.ContraABDRhombomere.cellID(All.ContraABDRhombomere.r5==1);...
                        All.ContraABDRhombomere.cellID(All.ContraABDRhombomere.r6==1);...
                        All.ContraABDRhombomere.cellID(All.ContraABDRhombomere.r7==1)];
                    
                    
All.ContraABDiOrdered = [All.ContraABDiRhombomere.cellID(All.ContraABDiRhombomere.r3==1);...
                        All.ContraABDiRhombomere.cellID(All.ContraABDiRhombomere.r4==1);...
                        All.ContraABDiRhombomere.cellID(All.ContraABDiRhombomere.r5==1);...
                        All.ContraABDiRhombomere.cellID(All.ContraABDiRhombomere.r6==1);...
                        All.ContraABDiRhombomere.cellID(All.ContraABDiRhombomere.r7==1)];

%%

clear MatIndex
clear MatOrder
clear connMat


%Integrator, Saccadic, vestibular, Abducens
  MatOrder = [All.IntegratorABD;All.IntegratorABDi;...
              All.SaccadicABD; All.SaccadicABDi;...
              All.VestibularABD;All.VestibularABDi;...
              All.ContraABD;All.ContraABDi;...
              All.RestABD;All.RestABDi;...
              [ABDr.cellID]';[ABDc.cellID]';[ABDIr.cellID]';[ABDIc.cellID]'];

 %MatOrder = [All.ContraABDOrdered;All.ContraABDiOrdered;...
 %             [ABDr.cellID]';[ABDc.cellID]';[ABDIr.cellID]';[ABDIc.cellID]'];
              
for i = 1:numel(MatOrder)
    MatIndex(i) = find(AllCells == MatOrder(i),1);
end

connMat = ConnMatrixPre(MatIndex,MatIndex);
%randomized = randi(385,1,385);
cspy(connMat,'Colormap',colorcet('L1','reverse',1),'Levels',255,'MarkerSize',15);
%cspy(connMat([randomized,386:end],[randomized,386:end]),'Colormap',colorcet('L1','reverse',1),'Levels',255,'MarkerSize',15);

%imagesc(connMat);
endABDInt = find(All.IntegratorABD(end) == AllCells(MatIndex));
endABDiInt = find(All.IntegratorABDi(end) == AllCells(MatIndex));
endSacABD = find(All.SaccadicABD(end) == AllCells(MatIndex));
endSACABDi = find(All.SaccadicABDi(end) == AllCells(MatIndex));
endVestABD = find(All.VestibularABD(end) == AllCells(MatIndex));
endVestABDi = find(All.VestibularABDi(end) == AllCells(MatIndex));
endContraABD = find(All.ContraABD(end) == AllCells(MatIndex));
endContraABDi = find(All.ContraABDi(end) == AllCells(MatIndex));
endRestABD = find(All.RestABD(end) == AllCells(MatIndex));
endRestABDi = find(All.RestABDi(end) == AllCells(MatIndex));

endABD = find(77296 == AllCells(MatIndex));
hold on;

line([0,size(connMat,1)],[endABDInt,endABDInt],'color','k');
line([0,size(connMat,1)],[endABDiInt,endABDiInt],'color','k');
line([0,size(connMat,1)],[endSacABD,endSacABD],'color','k');
line([0,size(connMat,1)],[endSACABDi(1),endSACABDi(1)],'color','k');
line([0,size(connMat,1)],[endVestABD,endVestABD],'color','k');
line([0,size(connMat,1)],[endVestABDi,endVestABDi],'color','k');
line([0,size(connMat,1)],[endContraABD,endContraABD],'color','k');
line([0,size(connMat,1)],[endContraABDi,endContraABDi],'color','k');
line([0,size(connMat,1)],[endRestABD,endRestABD],'color','k');
line([0,size(connMat,1)],[endRestABDi,endRestABDi],'color','k');
line([0,size(connMat,1)],[endABD,endABD],'color','k');



line([endABDInt,endABDInt],[0,size(connMat,1)],'color','k');
line([endABDiInt,endABDiInt],[0,size(connMat,1)],'color','k');
line([endSacABD,endSacABD],[0,size(connMat,1)],'color','k');
line([endSACABDi(1),endSACABDi(1)],[0,size(connMat,1)],'color','k');
line([endVestABD,endVestABD],[0,size(connMat,1)],'color','k');
line([endVestABDi,endVestABDi],[0,size(connMat,1)],'color','k');
line([endContraABD,endContraABD],[0,size(connMat,1)],'color','k');
line([endContraABDi,endContraABDi],[0,size(connMat,1)],'color','k');
line([endRestABD,endRestABD],[0,size(connMat,1)],'color','k');
line([endRestABDi,endRestABDi],[0,size(connMat,1)],'color','k');
line([endABD,endABD],[0,size(connMat,1)],'color','k');
