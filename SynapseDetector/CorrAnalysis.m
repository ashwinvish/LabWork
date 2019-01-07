
vars = {'corr_func','meanCorr','corr_funcDBX','meanCorrDBX','corr_funcALX','meanCorrALX','corr_funcBARHL','meanCorrBARHL'};
clear(vars{:});


t = -2:0.05:7;

windowWidth = 1; % in seconds
windows  = 1:windowWidth/0.05:181; % in frames

% All Cells
for i =2:size(windows,2)
    temp = corrcoef(Firing(windows(i-1):windows(i),:));
    temp = temp - diag(diag(temp));
    corr_func{i} = temp;
    meanCorr(i) = mean(tril(corr_func{i}(:)));
end

subplot(4,4,1)
shadedErrorBar(t,mean(Firing(:,OriginalCellOrderDBX),2),std(Firing(:,OriginalCellOrderDBX),[],2),{'color',colors(1,:)},0.5);
hold on
shadedErrorBar(t,mean(Firing(:,OriginalCellOrderALX),2),std(Firing(:,OriginalCellOrderALX),[],2),{'color',colors(2,:)},0.5);
shadedErrorBar(t,mean(Firing(:,OriginalCellOrderBARHL),2),std(Firing(:,OriginalCellOrderBARHL),[],2),{'color',colors(3,:)},0.5);

set(gca,'XLim',[-2,7]);
box off;

%%
%DBXCells

for i =2:size(windows,2)
    temp = corrcoef(Firing(windows(i-1):windows(i),OriginalCellOrderDBX));
    temp = temp - diag(diag(temp));
    corr_funcDBX{i} = temp;
    meanCorrDBX(i) = mean(tril(corr_funcDBX{i}(:)));
end

for i = 2:10
    corr_DBX1(:,i) = corr_funcDBX{i}(1,:)';
end

%%

% ALX

for i =2:size(windows,2)
    temp = corrcoef(Firing(windows(i-1):windows(i),OriginalCellOrderALX));
    temp = temp - diag(diag(temp));
    corr_funcALX{i} = temp;
    meanCorrALX(i) = mean(tril(corr_funcALX{i}(:)));
end

% Barhl

for i =2:size(windows,2)
    temp = corrcoef(Firing(windows(i-1):windows(i),OriginalCellOrderBARHL));
    temp = temp - diag(diag(temp));
    corr_funcBARHL{i} = temp;
    meanCorrBARHL(i) = mean(tril(corr_funcBARHL{i}(:)));
end

colors = cbrewer('qual', 'Set1', 10);
subplot(4,4,2)
plot(1:size(windows,2),meanCorrDBX,'Marker','o','MarkerSize',10,'Color',colors(1,:),'MarkerFaceColor',colors(1,:),'MarkerEdgeColor','none','LineWidth',2);
hold on;
plot(1:size(windows,2),meanCorrALX,'Marker','o','MarkerSize',10,'Color',colors(2,:),'MarkerFaceColor',colors(2,:),'MarkerEdgeColor','none','LineWidth',2);
plot(1:size(windows,2),meanCorrBARHL,'Marker','o','MarkerSize',10,'Color',colors(3,:),'MarkerFaceColor',colors(3,:),'MarkerEdgeColor','none','LineWidth',2);
box off;
xticks(1:2:size(windows,2));
xticklabels({'-1.950','0.05','2.05','4.05' ,'6.05'});
legend({'DBX','ALX','BARHL'});

