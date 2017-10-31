PrimaryCilia = [allTrees{1}{1,2}{1,3};
    allTrees{2}{1,2}{1,3};
    [5*12069, 5*20833, 45*70];
    allTrees{4}{1,2}{1,3};
    allTrees{5}{1,2}{1,3};
    allTrees{6}{1,3}{1,3};
    allTrees{7}{1,2}{1,3};
    allTrees{8}{1,8}{1,3};
    allTrees{9}{1,2}{1,3};
    allTrees{10}{1,2}{1,3};
    allTrees{11}{1,3}{1,3};
    allTrees{12}{1,3}{1,3};
    allTrees{13}{1,2}{1,3};
    allTrees{14}{1,3}{1,3};
    allTrees{15}{1,3}{1,3};
    allTrees{16}{1,2}{1,3};
    [5*19711, 5* 23075, 45*436];
    allTrees{18}{1,5}{1,3};
    allTrees{19}{1,3}{1,3};
    allTrees{20}{1,4}{1,3};
    allTrees{21}{1,3}{1,3};
    allTrees{22}{1,3}{1,3}];

AlxCiliaAngle = [];
DbxCiliaAngle = [];
TransCiliaAngle = [];
BarhlCiliaAngle = [];

%figure();
for i = 1:numel(cellIDs)
    if ismember(cellIDs{i}, cellIDsAlx) ==1
        %subplot(221)
        h1 = compass(PrimaryCilia(i,1)/1000-CellSoma(i,1)/1000, PrimaryCilia(i,2)/1000-CellSoma(i,2)/1000);
        temp =  CiliaAngle(PrimaryCilia(i,1)/1000-CellSoma(i,1)/1000, PrimaryCilia(i,2)/1000-CellSoma(i,2)/1000);
        h1.LineWidth = 2;
        h1.Color = calx;
        hold on;
        AlxCiliaAngle = [AlxCiliaAngle, temp];
        clear temp;
    end
end
 %axis([min(h1.XData) max(h1.XData), min(h1.YData), max(h1.YData)]);
 set(gca, 'FontName', 'Arial', 'FontSize', 40);
 
figure();
for i = 1:numel(cellIDs)
    if ismember(cellIDs{i}, cellIDsTrans) ==1
        %subplot(222)
        h2 = compass(PrimaryCilia(i,1)/1000-CellSoma(i,1)/1000, PrimaryCilia(i,2)/1000-CellSoma(i,2)/1000);
        temp = CiliaAngle(PrimaryCilia(i,1)/1000-CellSoma(i,1)/1000, PrimaryCilia(i,2)/1000-CellSoma(i,2)/1000);
        h2.LineWidth = 2;
        h2.Color = ctrans;
        hold on;
        TransCiliaAngle = [TransCiliaAngle,temp];
        clear temp;
    end
end
% axis([min(h2.XData) max(h2.XData), min(h2.YData), max(h2.YData)])
set(gca, 'FontName', 'Arial', 'FontSize', 40);

figure();
for i = 1:numel(cellIDs)

    if ismember(cellIDs{i}, cellIDsDbx) ==1
        %subplot(223)
        h3 = compass(PrimaryCilia(i,1)/1000-CellSoma(i,1)/1000, PrimaryCilia(i,2)/1000-CellSoma(i,2)/1000);
        temp = CiliaAngle(PrimaryCilia(i,1)/1000-CellSoma(i,1)/1000, PrimaryCilia(i,2)/1000-CellSoma(i,2)/1000);
        h3.LineWidth = 2;
        h3.Color = cdbx;
        hold on;
        DbxCiliaAngle = [DbxCiliaAngle, temp];
        clear temp;
    end
end
%axis([min(h3.XData) max(h3.XData), min(h3.YData), max(h3.YData)])
set(gca, 'FontName', 'Arial', 'FontSize', 40);

figure();
for i = 1:numel(cellIDs)
    if ismember(cellIDs{i}, cellIDsL) ==1
        %subplot(224)
        h4 = compass(PrimaryCilia(i,1)/1000-CellSoma(i,1)/1000, PrimaryCilia(i,2)/1000-CellSoma(i,2)/1000);
        temp = CiliaAngle(PrimaryCilia(i,1)/1000-CellSoma(i,1)/1000, PrimaryCilia(i,2)/1000-CellSoma(i,2)/1000);
        h4.LineWidth = 2;
        h4.Color = cbarhl;
        hold on;
        BarhlCiliaAngle = [BarhlCiliaAngle, temp];
        clear temp;
    end
 
end
  % axis([min(h4.XData) max(h4.XData), min(h4.YData), max(h4.YData)])
set(gca, 'FontName', 'Arial', 'FontSize', 40);

% get all angles

figure();

h = boxplot([AlxCiliaAngle'; TransCiliaAngle'; DbxCiliaAngle'; BarhlCiliaAngle'],...
    [ones(size(AlxCiliaAngle,2),1); 2*ones(size(TransCiliaAngle,2),1); 3*ones(size(DbxCiliaAngle,2),1); 4*ones(size(BarhlCiliaAngle,2),1)],'Notch','off', 'Symbol', 'ko',...
    'Colors',[calx;ctrans;cdbx;cbarhl],'OutlierSize', 25);
set(h,{'linew'},{4});
set(findobj(gcf,'LineStyle','--'),'LineStyle','-');
h1 = findobj(gca,'tag','Median');
set(h1,'Color','k');
set(gca,'FontName', 'Arial', 'FontSize', 40, 'LineWidth',4);
box off;


    
    
    
    





