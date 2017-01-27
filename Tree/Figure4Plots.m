% all alx cells
PlotClass(cellIDsAlx,allTrees,cellIDs, calx);
h1 = gcf;
set(h1,'color','none', 'units','normalized','outerposition',[0 0 1 1]);
axis vis3d;
[h2,h3] = PlotViews(h1);
close(figure(4));


export_fig(h1,sprintf('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/ZfishpaperFigures/NewFigure4/AllAlxXY.png'),'-r300','-transparent');
export_fig(h2,sprintf('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/ZfishpaperFigures/NewFigure4/AllAlxYZ.png'),'-r300','-transparent');
export_fig(h3,sprintf('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/ZfishpaperFigures/NewFigure4/AllAlxXZ.png'),'-r300','-transparent');

pause(100);
close all;
clear col;

% all Trans cells
PlotClass(cellIDsTrans,allTrees,cellIDs, ctrans);
h1 = gcf;
set(h1,'color','none', 'units','normalized','outerposition',[0 0 1 1]);
axis vis3d;
[h2,h3] = PlotViews(h1);
close(figure(4));


export_fig(h1,sprintf('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/ZfishpaperFigures/NewFigure4/AllTransXY.png'),'-r300','-transparent');
export_fig(h2,sprintf('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/ZfishpaperFigures/NewFigure4/AllTransYZ.png'),'-r300','-transparent');
export_fig(h3,sprintf('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/ZfishpaperFigures/NewFigure4/AllTransXZ.png'),'-r300','-transparent');

pause(100);
close all;
clear col;

% all Dbx cells
PlotClass(cellIDsDbx,allTrees,cellIDs, cdbx);
h1 = gcf;
set(h1,'color','none', 'units','normalized','outerposition',[0 0 1 1]);
axis vis3d;
[h2,h3] = PlotViews(h1);
close(figure(4));

export_fig(h1,sprintf('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/ZfishpaperFigures/NewFigure4/AllDbxXY.png'),'-r300','-transparent');
export_fig(h2,sprintf('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/ZfishpaperFigures/NewFigure4/AllDbxYZ.png'),'-r300','-transparent');
export_fig(h3,sprintf('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/ZfishpaperFigures/NewFigure4/AllDbxXZ.png'),'-r300','-transparent');

pause(100);
close all;
clear col;

% all barhl cells

PlotClass(cellIDsL,allTrees,cellIDs, cbarhl);
h1 = gcf;
set(h1,'color','none', 'units','normalized','outerposition',[0 0 1 1]);
axis vis3d;
[h2,h3] = PlotViews(h1);
close(figure(4));


export_fig(h1,sprintf('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/ZfishpaperFigures/NewFigure4/AllBarhlXY.png'),'-r300','-transparent');
export_fig(h2,sprintf('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/ZfishpaperFigures/NewFigure4/AllBarhlYZ.png'),'-r300','-transparent');
export_fig(h3,sprintf('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/ZfishpaperFigures/NewFigure4/AllBarhlXZ.png'),'-r300','-transparent');

pause(100);
close all;
clear col;

