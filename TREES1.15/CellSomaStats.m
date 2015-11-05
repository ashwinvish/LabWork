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
    plot(MeanRCEucDist(i)./1000,tempdiffRC,'o', 'MarkerEdgeColor', [0.7,0.7,0.7], 'MarkerFaceColor',[0.7,0.7,0.7], 'MarkerSize', 20 );
    hold on;
end
figure(1);
plot(MeanRCEucDist./1000,RhoDiffRC,'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor','none', 'MarkerSize', 20 );
plot([MeanRCEucDist./1000;MeanRCEucDist./1000], [RhoDiffRC-RhoDiffRC_SD; RhoDiffRC+RhoDiffRC_SD], 'Color','k','LineWidth',2);
xlabel('Pairwise RC Eucledian distance in \mum', 'FontName', 'Arial', 'FontSize', 20);
ylabel('Pairwise difference in persistence measure \rho', 'FontName', 'Arial', 'FontSize', 20);
set(gca, 'FontName', 'Arial', 'FontSize', 20);
set(gcf,'color','w');
box off;


X = [ones(length(MeanRCEucDist./1000),1) MeanRCEucDist'./1000];
y =RhoDiffRC';
b = X\y;
yCalc2 = X*b;
Rsq2 = 1 - sum((y - yCalc2).^2)/sum((y - mean(y)).^2);
plot(MeanRCEucDist'./1000,yCalc2,'-r','LineWidth',2);
text(max(MeanRCEucDist'./1000),max(yCalc2), sprintf('R^2 = %0.2f',Rsq2), 'FontName', 'Arial', 'FontSize', 20 );
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
    plot(MeanDVEucDist(i)./1000,tempdiffDV,'o', 'MarkerEdgeColor', [0.7,0.7,0.7], 'MarkerFaceColor',[0.7,0.7,0.7], 'MarkerSize', 20   );
    hold on;
end
figure (2);
plot(MeanDVEucDist./1000,RhoDiffDV,'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor','none', 'MarkerSize', 20 );
set(gca,'XLim',[0 23]);
plot([MeanDVEucDist./1000 ; MeanDVEucDist./1000], [RhoDiffDV-RhoDiffDV_SD; RhoDiffDV+RhoDiffDV_SD], 'Color','k','LineWidth',2);
xlabel('Pairwise DV Eucledian distance in \mum', 'FontName', 'Arial', 'FontSize', 20);
ylabel('Pairwise difference in persistence measure \rho', 'FontName', 'Arial', 'FontSize', 20);
box off

X = [ones(length(MeanDVEucDist./1000),1) MeanDVEucDist'./1000];
y =RhoDiffDV';
b = X\y;
yCalc2 = X*b;
Rsq2 = 1 - sum((y - yCalc2).^2)/sum((y - mean(y)).^2);
plot(MeanDVEucDist'./1000,yCalc2,'-r','LineWidth',2);
text(max(MeanDVEucDist'./1000),max(yCalc2), sprintf('R^2 = %0.2f',Rsq2) , 'FontName', 'Arial', 'FontSize', 20);
set(gca, 'FontName', 'Arial', 'FontSize', 20);
set(gcf,'color','w');
axis square;
clear temp;

%%

%Distribution of Persistence times along axis
calx = [1 0.5 0.3];
cdbx = [1 0.3 1];
cbarhl = [0.3 0.5 1];
RhoAlx = [];
RhoDbx =[];
RhoBarhl = [];
[y,I] = sort(rho);


for i = 1:length(cellIDs)
    if ismember(cellIDs(I(i)),cellIDsAlx) == 1
        BarCMap = calx;
        RhoAlx = [RhoAlx,rho(I(i))];
    elseif ismember(cellIDs(I(i)),cellIDsDbx) == 1
        BarCMap = cdbx;
        RhoDbx = [RhoDbx,rho(I(i))];
    else
        BarCMap= cbarhl;
        RhoBarhl = [RhoBarhl,rho(I(i))];
    end
    
    h = plot(i,rho(I(i)),'o','MarkerSize',20, 'MarkerFaceColor',BarCMap, 'MarkerEdgeColor','none');
    %set(h, 'FaceColor', BarCMap);
    hold on;
    
end

xlabel('Neuron #', 'FontName', 'Arial', 'FontSize', 20);
ylabel('Persistence time Measure \rho', 'FontName', 'Arial', 'FontSize', 20);
set(gca, 'FontName', 'Arial', 'FontSize', 20);
set(gcf,'color','w');
axis square;
box off


figure();
MeanRhos = [mean(RhoAlx); mean(RhoDbx); mean(RhoBarhl)];
Colors = [calx;cdbx;cbarhl];
 for i = 1:3 
    h =  bar(i,MeanRhos(i));
    set(h,'FaceColor',Colors(i,:));
    hold on;
 end
 box off
 axis square
 
 set(gca,'XTick', [1:3],'XTickLabel', {'group1'; 'group2'; 'group3'}, 'FontName', 'Arial', 'FontSize', 20);
 set(gca, 'FontName', 'Arial', 'FontSize', 20);
 set(gcf,'color','w');
 ylabel('Persistence time Measure \rho', 'FontName', 'Arial', 'FontSize', 20);
 


