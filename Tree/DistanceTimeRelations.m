RCSoma = [CellSoma(:,1),CellSoma(:,2)]; % considering only the ML,RC coordinates
RCSomapDist = tril(squareform(pdist(RCSoma)),-1);
%cmap = parula(44);
index =1;
pairsDist = [];
pairsRhoDiff = [];

% RC-ML
figure();
for i = 1:numel(cellIDs)
        cmap = parula(22); 
subplot(3,8,i);
    for j= 1:numel(cellIDs)
        pairsDist = [pairsDist;(RCSoma(1,:)-RCSoma(j,:))];
        pairsRhoDiff = [pairsRhoDiff;abs(rho(i) - rho(j))];
        %scatter(pairsDist(1),pairsDist(2),20,parula(length(pairsRhoDiff)));
        %hold on;
        plot(pairsDist(:,1)/1000,pairsDist(:,2)/1000,'o');
        line(pairsDist(:,1)/1000,pairsDist(:,2)/1000,'Color',cmap(j,:));
        hold on;
    end
    
    cmap = parula(length(pairsRhoDiff));
    %scatter(pairsDist(:,1)/1000,pairsDist(:,2)/1000,70,pairsRhoDiff,'filled');
    %plot(pairsDist(:,1)/1000,pairsDist(:,2)/1000,'-o', 'Color', cmap(i,:));
    set(gca,'XLim',[-100,100],'YLim',[-100,100],'YDir','reverse');
    str = sprintf('%s', cellIDs{i});
    title(str);
    colormap jet;
    hold on;
    PlotAxisAtOrigin(0,0);
    axis square;
    pairsDist = [];
    pairsRhoDiff = [];
    %hold on;
end

%RCAxis only
figure();
pairsDist = [];
pairsRhoDiff = [];
for i = 1:numel(cellIDs)
    for j= 1:numel(cellIDs)
        pairsDist = [pairsDist;(RCSoma(i,2)-RCSoma(j,2))];
        pairsRhoDiff = [pairsRhoDiff;abs(rho(i) - rho(j))];
        %scatter(pairsDist(1),pairsDist(2),20,parula(length(pairsRhoDiff)));
        %hold on;
    end
    subplot(3,8,i);
    scatter(zeros(length(pairsDist),1),pairsDist/1000,70,pairsRhoDiff,'filled');
    set(gca,'XLim',[-1,1],'YLim',[-100,100],'YDir','reverse');
    str = sprintf('%s', cellIDs{i});
    title(str);
    colormap jet;
    hold on;
    PlotAxisAtOrigin(0,0);
    %axis square;
    pairsDist = [];
    pairsRhoDiff = [];
    %hold on;
end

%%DV Axis only

figure();
DVSoma = CellSoma(:,3);
pairsDist = [];
pairsRhoDiff = [];
for i = 1:numel(cellIDs)
    for j= 1:numel(cellIDs)
        pairsDist = [pairsDist;(DVSoma(i)-DVSoma(j))];
        pairsRhoDiff = [pairsRhoDiff;abs(rho(i) - rho(j))];
        %scatter(pairsDist(1),pairsDist(2),20,parula(length(pairsRhoDiff)));
        %hold on;
    end
    subplot(3,8,i);
    scatter(zeros(length(pairsDist),1),pairsDist/1000,70,pairsRhoDiff,'filled');
    set(gca,'XLim',[-1,1],'YLim',[-60,60],'YDir','reverse');
    hold on;
    str = sprintf('%s', cellIDs{i});
    title(str);
    colormap jet;
    PlotAxisAtOrigin(0,0);
    %axis square;
    pairsDist = [];
    pairsRhoDiff = [];
    %hold on;
end


