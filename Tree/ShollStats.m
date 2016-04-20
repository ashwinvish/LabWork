
radius = 10000;

for i = 1:size(cellIDs,2)
    tree{i} = load_tree([cellIDs{i} , '_WithTags.swc']);
    shollIntersections{i}= sholl_tree(tree{i},radius); % concentric circles every radius nm
    plot(1:radius:radius*length(shollIntersections{i}), shollIntersections{i});
    hold on;
end
radius = radius/1000; % convert to um
hold off;
%[SortLength, I] = sort(cell2mat(allRawLength(:))/1000);
[SortShollIntersetctions, I] = sort(cellfun(@sum, shollIntersections));


ShollAlx = [];
ShollTrans = [];
ShollDbx = [];
ShollBarhl = [];

indexAlx = 1;
indexTrans = 1;
indexDbx = 1;
indexBarhl = 1;

figure()


for i = 1:length(cellIDs)
    if ismember(cellIDs(I(i)),cellIDsAlx) == 1
        BarCMap = calx;
        ShollAlx = [ShollAlx,SortShollIntersetctions(i)];
        ShollAlxIntersections(indexAlx,1:length(shollIntersections{I(i)})) = shollIntersections{I(i)};
        indexAlx = indexAlx+1;
        
    elseif ismember(cellIDs(I(i)),cellIDsTrans) == 1
        BarCMap = ctrans;
        ShollTrans = [ShollTrans,SortShollIntersetctions(i)];
        ShollTransIntersections(indexTrans,1:length(shollIntersections{I(i)})) = shollIntersections{I(i)};
        indexTrans= indexTrans+1;
        
    elseif ismember(cellIDs(I(i)),cellIDsDbx) == 1
        BarCMap = cdbx;
        ShollDbx = [ShollDbx,SortShollIntersetctions(i)];
        ShollDbxIntersections(indexDbx,1:length(shollIntersections{I(i)})) = shollIntersections{I(i)};
        indexDbx= indexDbx+1;
        
    else
        BarCMap= cbarhl;
        ShollBarhl = [ShollBarhl,SortShollIntersetctions(i)];
        ShollBarhlIntersections(indexBarhl,1:length(shollIntersections{I(i)})) = shollIntersections{I(i)};
        indexBarhl = indexBarhl+1;
    end
    h = bar(i,SortShollIntersetctions(i),'FaceColor', BarCMap);
    hold on;
    
end

figure();

%     %plot(1:radius:radius*length(ShollAlxIntersections(i,:)), ShollAlxIntersections(i,:),'o', 'MarkerFaceColor', calx );
%     plot(1:radius:radius*length(mean(ShollAlxIntersections)), mean(ShollAlxIntersections),'color', calx, 'LineWidth', 2);
%     hold on;
%     plot([1:radius:radius*length(mean(ShollAlxIntersections));1:radius:radius*length(mean(ShollAlxIntersections))],[mean(ShollAlxIntersections)-std(ShollAlxIntersections);mean(ShollAlxIntersections)+std(ShollAlxIntersections)],...
%         'color', calx, 'LineWidth', 2);
%     %plot(1:radius:radius*length(ShollDbxIntersections(i,:)), ShollDbxIntersections(i,:),'o', 'MarkerFaceColor', cdbx);
%     plot(1:radius:radius*length(mean(ShollDbxIntersections)), mean(ShollDbxIntersections),'color', cdbx, 'LineWidth', 2);
%     %plot(1:radius:radius*length(ShollBarhlIntersections(i,:)), ShollBarhlIntersections(i,:),'o', 'MarkerFaceColor', cbarhl);
%     plot(1:radius:radius*length(mean(ShollBarhlIntersections)), mean(ShollBarhlIntersections),'color', cbarhl, 'LineWidth', 2);


shadedErrorBar((1:radius:radius*length(mean(ShollAlxIntersections))),ShollAlxIntersections,{@mean,@std},{'-','color',calx,'markerfacecolor',calx,'LineWidth',4},1);  
hold on;
shadedErrorBar((1:radius:radius*length(mean(ShollTransIntersections))),ShollTransIntersections,{@mean,@std},{'-','color',ctrans,'markerfacecolor',ctrans,'LineWidth',4},1);  
shadedErrorBar((1:radius:radius*length(mean(ShollDbxIntersections))),ShollDbxIntersections,{@mean,@std},{'-','color',cdbx,'markerfacecolor',cdbx,'LineWidth',4},1);  
shadedErrorBar((1:radius:radius*length(mean(ShollBarhlIntersections))),ShollBarhlIntersections,{@mean,@std},{'-','color',cbarhl,'markerfacecolor',cbarhl,'LineWidth',4},1);  


set(gca, 'FontSize',40, 'FontName', 'Arial', 'LineWidth',2);
xlabel('Length (\mum)', 'FontSize',40, 'FontName', 'Arial');
ylabel('Mean Sholl intersetcions', 'FontSize',40, 'FontName', 'Arial'); 
set(gca,'XLim', [0 220]);
set(gcf,'color','w');
axis square
box off;


