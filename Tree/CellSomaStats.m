%% RC axis only
RCSoma = [CellSoma(:,1),CellSoma(:,2)]; % considering only the x,y coordinates
RCpdist = tril(squareform(pdist(RCSoma)),-1);
clear RhoDiffRC;
clear RhoDiffRC_SD;
clear MeanRCEucDist;

index = 1;
steps = 1:10000:max(RCpdist(:));
for i = 1:length(steps)-1
    tempRC = find(RCpdist>steps(i) & RCpdist<steps(i+1));
    [m,n] = ind2sub(size(RCpdist),tempRC);
    tempdiffRC = abs(rho(m)-rho(n));
    RhoDiffRC(index) = mean(abs(rho(m)-rho(n)));
    RhoDiffRC_SD(index) = std(abs(rho(m)-rho(n)));
    MeanRCEucDist(index) = mean(RCpdist(tempRC));
    index = index+1;
    figure(1);
    plot(MeanRCEucDist(i)./1000,tempdiffRC,'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor',[0.9,0.9,0.9], 'MarkerSize', 25 );
    hold on;
end
figure(1);
plot(MeanRCEucDist./1000,RhoDiffRC,'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor','none', 'MarkerSize', 35 );
plot([MeanRCEucDist./1000;MeanRCEucDist./1000], [RhoDiffRC-RhoDiffRC_SD; RhoDiffRC+RhoDiffRC_SD], 'Color','k','LineWidth',2);
xlabel('Pairwise RC distance in \mum', 'FontName', 'Arial', 'FontSize', 40);
ylabel('Pairwise persistence difference \rho', 'FontName', 'Arial', 'FontSize', 40);
set(gca, 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
set(gcf,'color','w');
box off;


X = [ones(length(MeanRCEucDist./1000),1) MeanRCEucDist'./1000];
y =RhoDiffRC';
b = X\y;
yCalc2 = X*b;
Rsq2 = 1 - sum((y - yCalc2).^2)/sum((y - mean(y)).^2);
plot(MeanRCEucDist'./1000,yCalc2,'-r','LineWidth',2);
%text(max(MeanRCEucDist'./1000),max(yCalc2), sprintf('R^2 = %0.2f',Rsq2), 'FontName', 'Arial', 'FontSize', 20 );
[PearsonsCoeffRC, PvalRC] = corr(MeanRCEucDist',RhoDiffRC');
axis square;
clear temp;



% DV axis
clear MeanDVEucDist;
clear RhoDiffDV;
clear RhoDiffDV_SD;

DVpdist = tril(squareform(pdist(CellSoma(:,3))),-1);

index = 1;
steps = 1:2500:max(DVpdist(:));
for i = 1:length(steps)-1
    tempDV = find(DVpdist>steps(i) & DVpdist<steps(i+1));
    [m,n] = ind2sub(size(DVpdist),tempDV);
    tempdiffDV = abs(rho(m)-rho(n));    
    RhoDiffDV(index) = mean(abs(rho(m)-rho(n)));
    RhoDiffDV_SD(index) = std(abs(rho(m)-rho(n)));
    MeanDVEucDist(index) = mean(DVpdist(tempDV));
    index = index+1;
    figure(2);
    plot(MeanDVEucDist(i)./1000,tempdiffDV,'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor',[0.7,0.7,0.7], 'MarkerSize', 25   );
    hold on;
end
figure (2);
plot(MeanDVEucDist./1000,RhoDiffDV,'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor','none', 'MarkerSize', 35);
set(gca,'XLim',[0 23]);
plot([MeanDVEucDist./1000 ; MeanDVEucDist./1000], [RhoDiffDV-RhoDiffDV_SD; RhoDiffDV+RhoDiffDV_SD], 'Color','k','LineWidth',2);
xlabel('Pairwise DV distance in \mum', 'FontName', 'Arial', 'FontSize', 40);
ylabel('Pairwise persistence difference', 'FontName', 'Arial', 'FontSize', 40);
box off

X = [ones(length(MeanDVEucDist./1000),1) MeanDVEucDist'./1000];
y =RhoDiffDV';
b = X\y;
yCalc2 = X*b;
Rsq2 = 1 - sum((y - yCalc2).^2)/sum((y - mean(y)).^2);
plot(MeanDVEucDist'./1000,yCalc2,'-r','LineWidth',2);
%text(max(MeanDVEucDist'./1000),max(yCalc2), sprintf('R^2 = %0.2f',Rsq2) , 'FontName', 'Arial', 'FontSize', 20);
[PearsonCoeffDV, PvalDV] = corr(MeanDVEucDist',RhoDiffDV');
set(gca, 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
set(gcf,'color','w');
axis square;
clear temp;

%%
clear m;
clear n;
clear RCSoma;
clear RCpdist;
clear tempRC;
clear tempdiffRC;

RCSoma = [CellSoma(:,1),CellSoma(:,2)]; % considering only the x,y coordinates
RCpdist = tril(squareform(pdist(RCSoma)),-1);

% tempRC = find(RCpdist>1 & RCpdist< max(RCpdist(:)));

tempRC = find(RCpdist>1 & RCpdist<20000);

[m,n] = ind2sub(size(RCpdist),tempRC);
tempdiffRC = abs(rho(m)-rho(n));
corrcoef(RCpdist(tempRC)/1000,tempdiffRC)
plot(RCpdist(tempRC)/1000,tempdiffRC,'o','MarkerSize',25, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.7,0.7,0.7]);
box off;
set(gca, 'FontName', 'Arial', 'FontSize', 40, 'LineWidth',2);
set(gcf,'color','w');
axis square;

%%

%Distribution of Persistence times along axis
calx = [0.9655    0.6207    0.8621];%[1,0.5,0];
cdbx = [1.0000    0.7586    0.5172];%[1, 0, 1];
ctrans = [ 0.5517    0.6552    0.4828];%[0,1,0.5];
cbarhl = [0.6207    0.7586    1.0000];%[0, 0.5, 1];

RhoAlx = [];
RhoTrans= []
RhoDbx =[];
RhoBarhl = [];
[y,I] = sort(rho);


for i = 1:length(cellIDs)
    if ismember(cellIDs(I(i)),cellIDsAlx) == 1
        BarCMap(i,:) = calx;
        RhoAlx = [RhoAlx,rho(I(i))];
    elseif ismember(cellIDs(I(i)),cellIDsTrans) == 1
        BarCMap(i,:) = ctrans;
        RhoTrans = [RhoTrans,rho(I(i))];
        
    elseif ismember(cellIDs(I(i)),cellIDsDbx) == 1
        BarCMap(i,:) = cdbx;
        RhoDbx = [RhoDbx,rho(I(i))];
    else
        BarCMap(i,:) = cbarhl;
        RhoBarhl = [RhoBarhl,rho(I(i))];
    end
    
    h = plot(i,rho(I(i)),'o','MarkerSize',25, 'MarkerFaceColor',BarCMap(i,:), 'MarkerEdgeColor','k');
    %set(h, 'FaceColor', BarCMap);
    hold on;
    
end

xlabel('Neuron #', 'FontName', 'Arial', 'FontSize', 40);
ylabel('Persistence time Measure \rho', 'FontName', 'Arial', 'FontSize', 40);
set(gca, 'FontName', 'Arial', 'FontSize', 40, 'LineWidth', 2);
set(gcf,'color','w');
axis square;
box off


figure();
MeanRhos = [mean(RhoAlx); mean(RhoTrans); mean(RhoDbx); mean(RhoBarhl)];
Colors = [calx;ctrans;cdbx;cbarhl];
 for i = 1:4
    h =  bar(i,MeanRhos(i));
    set(h,'FaceColor',Colors(i,:));
    hold on;
 end
 box off
 axis square
 
 set(gca,'XTick', [1:4],'XTickLabel', {'Ipsi'; 'Ipsi-Contra';'Contra'; 'Unknown'}, 'FontName', 'Arial', 'FontSize', 20);
 set(gca, 'FontName', 'Arial', 'FontSize', 20);
 set(gcf,'color','w');
 ylabel('Persistence time Measure \rho', 'FontName', 'Arial', 'FontSize', 20);
 


