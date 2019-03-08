%% Make plot for axon organization


load('AllCells.mat');

%

ContraLogic = isContra(AllCells);
ConfirmedContra = AllCells(ContraLogic);
RemainingCells = setdiff(AllCells,ConfirmedContra);
ConfirmedRemainingCells = RemainingCells(isExistReRoot(RemainingCells));
ContraRhombomere = isRhombomere(ConfirmedContra);
IpsiRhombomere  = isRhombomere(ConfirmedRemainingCells);
save('ContraRhombomere.mat','ContraRhombomere');
save('IpsiRhombomere.mat','IpsiRhombomere');


%% make nice plots

% contra cells only plot based on RC position

Colors = cmap('R2','N',5);

transform_swc_AV(ContraRhombomere.cellID(find(ContraRhombomere.r3 ==1)),Colors(1,:),[],true,false);
hold on;
transform_swc_AV(ContraRhombomere.cellID(find(ContraRhombomere.r4 ==1)),Colors(2,:),[],false,false);
transform_swc_AV(ContraRhombomere.cellID(find(ContraRhombomere.r5 ==1)),Colors(3,:),[],false,false);
transform_swc_AV(ContraRhombomere.cellID(find(ContraRhombomere.r6 ==1)),Colors(4,:),[],false,false);
transform_swc_AV(ContraRhombomere.cellID(find(ContraRhombomere.r7 ==1)),Colors(5,:),[],false,false);

colormap(Colors)

colorbar('Ticks',[0.1,0.3,0.5,0.7,0.9],...
         'TickLabels',{'r3','r4','r5','r6','r7'},'Location','east','Position',[.35 .5 ,0.025 .1]);
     
export_fig('/Users/ashwin/Desktop/RhombomereContraOrganization.png','-r300','-transparent');

%% RC organization by location
[id,order] = sort(ContraRhombomere.rootNode(:,2));
Colors = cmap('R2','N',size(order,1));

transform_swc_AV(ContraRhombomere.cellID(order),Colors,[],true,false);


%% DV by rhombomere

[id,order] = sort(ContraRhombomere.rootNode(:,3));
Colors = cmap('R2','N',size(order,1));

transform_swc_AV(ContraRhombomere.cellID(order),Colors,[],true,false);

%% organize by rhombomere

R3 = ContraRhombomere.cellID(ContraRhombomere.r3 ==1);
R3index = find(ContraRhombomere.r3 ==1);
[id,orderRC] = sort(ContraRhombomere.rootNode(R3index,2));
[id,orderDV] = sort(ContraRhombomere.rootNode(R3index,3));

Colors = cmap('R2','N',size(orderRC,1));

figure('units','normalized','outerposition',[0 0 1 1]);
transform_swc_AV(ContraRhombomere.cellID(R3index(orderRC)),Colors,[],true,false);
export_fig('/Users/ashwin/Desktop/RhombomereContraOrganization_RC_R3.png','-r300','-transparent');
close all;

figure('units','normalized','outerposition',[0 0 1 1]);
transform_swc_AV(ContraRhombomere.cellID(R3index(orderDV)),Colors,[],false,false);
view (90,0);
export_fig('/Users/ashwin/Desktop/RhombomereContraOrganization_DV_R3.png','-r300','-transparent');


% R4

R4 = ContraRhombomere.cellID(ContraRhombomere.r4 ==1);
R4index = find(ContraRhombomere.r4 ==1);
[id,orderRC] = sort(ContraRhombomere.rootNode(R4index,2));
[id,orderDV] = sort(ContraRhombomere.rootNode(R4index,3));

Colors = cmap('R2','N',size(orderRC,1));

figure('units','normalized','outerposition',[0 0 1 1]);
transform_swc_AV(ContraRhombomere.cellID(R4index(orderRC)),Colors,[],true,false);
export_fig('/Users/ashwin/Desktop/RhombomereContraOrganization_RC_R4.png','-r300','-transparent');
close all

figure('units','normalized','outerposition',[0 0 1 1]);
transform_swc_AV(ContraRhombomere.cellID(R4index(orderDV)),Colors,[],false,false);
view(90,0);
export_fig('/Users/ashwin/Desktop/RhombomereContraOrganization_DV_R4.png','-r300','-transparent');
close all;

% R5

R5 = ContraRhombomere.cellID(ContraRhombomere.r5 ==1);
R5index = find(ContraRhombomere.r5 ==1);
[id,orderRC] = sort(ContraRhombomere.rootNode(R5index,2));
[id,orderDV] = sort(ContraRhombomere.rootNode(R5index,3));

Colors = cmap('R2','N',size(orderRC,1));

figure('units','normalized','outerposition',[0 0 1 1]);
transform_swc_AV(ContraRhombomere.cellID(R5index(orderRC)),Colors,[],true,false);
export_fig('/Users/ashwin/Desktop/RhombomereContraOrganization_RC_R5.png','-r300','-transparent');
close all;

figure('units','normalized','outerposition',[0 0 1 1]);
transform_swc_AV(ContraRhombomere.cellID(R5index(orderDV)),Colors,[],false,false);
view (90,0);
export_fig('/Users/ashwin/Desktop/RhombomereContraOrganization_DV_R5.png','-r300','-transparent');
close all;

%R6

R6 = ContraRhombomere.cellID(ContraRhombomere.r6 ==1);
R6index = find(ContraRhombomere.r6 ==1);
[id,orderRC] = sort(ContraRhombomere.rootNode(R6index,2));
[id,orderDV] = sort(ContraRhombomere.rootNode(R6index,3));

Colors = cmap('R2','N',size(orderRC,1));

figure('units','normalized','outerposition',[0 0 1 1]);
transform_swc_AV(ContraRhombomere.cellID(R6index(orderRC)),Colors,[],true,false);
export_fig('/Users/ashwin/Desktop/RhombomereContraOrganization_RC_R6.png','-r300','-transparent');
close all;

figure('units','normalized','outerposition',[0 0 1 1]);
transform_swc_AV(ContraRhombomere.cellID(R6index(orderDV)),Colors,[],false,false);
view(90,0);
export_fig('/Users/ashwin/Desktop/RhombomereContraOrganization_DV_R6.png','-r300','-transparent');
close all;

% R7

R7 = ContraRhombomere.cellID(ContraRhombomere.r7 ==1);
R7index = find(ContraRhombomere.r7 ==1);
[id,orderRC] = sort(ContraRhombomere.rootNode(R7index,2));
[id,orderDV] = sort(ContraRhombomere.rootNode(R7index,3));

Colors = cmap('R2','N',size(orderRC,1));

figure('units','normalized','outerposition',[0 0 1 1]);
transform_swc_AV(ContraRhombomere.cellID(R7index(orderRC)),Colors,[],true,false);
export_fig('/Users/ashwin/Desktop/RhombomereContraOrganization_RC_R7.png','-r300','-transparent');
close all;

figure('units','normalized','outerposition',[0 0 1 1]);
transform_swc_AV(ContraRhombomere.cellID(R7index(orderDV)),Colors,[],false,false);
view(90,0)
export_fig('/Users/ashwin/Desktop/RhombomereContraOrganization_DV_R7.png','-r300','-transparent');
close all;


%%  Ipsi Cells

%r3
R3ipsi = IpsiRhombomere.cellID(IpsiRhombomere.r3 ==1);
R4ipsi = IpsiRhombomere.cellID(IpsiRhombomere.r4 ==1);
R5ipsi = IpsiRhombomere.cellID(IpsiRhombomere.r5 ==1);
R6ipsi = IpsiRhombomere.cellID(IpsiRhombomere.r6 ==1);
R7ipsi = IpsiRhombomere.cellID(IpsiRhombomere.r7 ==1);

Colors = cmap('R2','N',5);

figure('units','normalized','outerposition',[0 0 1 1]);
transform_swc_AV(R3ipsi,Colors(1,:),[],true,false);
export_fig('/Users/ashwin/Desktop/RhombomereIpsiOrganization_R3.png','-r300','-transparent');
close all;

%r4

figure('units','normalized','outerposition',[0 0 1 1]);
transform_swc_AV(R4ipsi,Colors(2,:),[],true,false);
export_fig('/Users/ashwin/Desktop/RhombomereIpsiOrganization_R4.png','-r300','-transparent');
close all;

%r5

figure('units','normalized','outerposition',[0 0 1 1]);
transform_swc_AV(R5ipsi,Colors(3,:),[],true,false);
export_fig('/Users/ashwin/Desktop/RhombomereIpsiOrganization_R5.png','-r300','-transparent');
close all;

%r6

figure('units','normalized','outerposition',[0 0 1 1]);
transform_swc_AV(R6ipsi,Colors(4,:),[],true,false);
export_fig('/Users/ashwin/Desktop/RhombomereIpsiOrganization_R6.png','-r300','-transparent');
close all;

%r7


figure('units','normalized','outerposition',[0 0 1 1]);
transform_swc_AV(R7ipsi,Colors(5,:),[],true,false);
export_fig('/Users/ashwin/Desktop/RhombomereIpsiOrganization_R7.png','-r300','-transparent');
close all;

% all together

figure('units','normalized','outerposition',[0 0 1 1]);

transform_swc_AV(R3ipsi,Colors(1,:),[],true,false);
transform_swc_AV(R4ipsi,Colors(2,:),[],false,false);
transform_swc_AV(R5ipsi,Colors(3,:),[],false,false);
transform_swc_AV(R6ipsi,Colors(4,:),[],false,false);
transform_swc_AV(R7ipsi,Colors(5,:),[],false,false);
export_fig('/Users/ashwin/Desktop/RhombomereIpsiOrganization_ALL.png','-r300','-transparent');
close all;


