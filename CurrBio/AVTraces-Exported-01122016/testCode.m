for i = 1:length(Small)
    if ismember(i,SmallAlx) ==1
        [a,h1,h2] = plotyy(1:386,peaksXY',Small(i,1),Small(i,3),'line','scatter');
        set(h1,'LineStyle','none');
        set(h2,'MarkerFaceColor',calx,'MarkerEdgeColor','none');
        set(a(2), 'XDir','reverse', 'YDir', 'reverse','YLim',[0,550], 'YTick',[], 'XTick',[]);
        hold on;
    elseif ismember(i, SmallTrans) ==1
        [a,h1,h2] = plotyy(1:386,peaksXY',Small(i,1),Small(i,3),'line','scatter');
        set(h1,'LineStyle','none');
        set(h2,'MarkerFaceColor',ctrans,'MarkerEdgeColor','none');
        set(a(2), 'XDir','reverse', 'YDir', 'reverse','YLim',[0,550], 'YTick',[], 'XTick',[]);
         hold on;
    elseif ismember( i, SmallDbx) ==1
        [a,h1,h2] = plotyy(1:386,peaksXY',Small(i,1),Small(i,3),'line','scatter');
        set(h1,'LineStyle','none');
        set(h2,'MarkerFaceColor',cdbx,'MarkerEdgeColor','none');
        set(a(2), 'XDir','reverse', 'YDir', 'reverse','YLim',[0,550], 'YTick',[], 'XTick',[]);
         hold on;
    else
        [a,h1,h2] = plotyy(1:386,peaksXY',Small(i,1),Small(i,3),'line','scatter')
        set(h1,'LineStyle','none');
        set(h2,'MarkerFaceColor',cbarhl,'MarkerEdgeColor','none');
        set(a(2), 'XDir','reverse', 'YDir', 'reverse','YLim',[0,550], 'YTick',[], 'XTick',[]);
         hold on;
    end
end
[a,h1,h2] = plotyy(1:386,peaksXY',0,0,'line','scatter');
set(h1,'LineStyle','-');
set(h2,'MarkerFaceColor','none','MarkerEdgeColor','none');
set(a(1),'XDir','reverse', 'YDir', 'reverse');
set(a(2), 'XDir','reverse', 'YDir', 'reverse','YLim',[0,550], 'YTick',[], 'XTick',[]);