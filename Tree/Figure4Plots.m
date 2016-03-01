% all alx cells
col = distinguishable_colors(numel(cellIDsAlx),calx);
DisplayTree(allTrees{4},[1],false, [eval([cellIDs{4},'_axon'])],col(1,:),allPreSynapse{4}, allPostSynapse{4});
DisplayTree(allTrees{6},[1],false, [eval([cellIDs{6},'_axon'])],col(3,:),allPreSynapse{6}, allPostSynapse{6});
set(gca,'color',[calx, 0.2]);
set(gcf,'color','none', 'units','normalized','outerposition',[0 0 1 1]);
h1 = gcf;
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
col = distinguishable_colors(numel(cellIDsTrans),ctrans);
DisplayTree(allTrees{7},[1],false, [eval([cellIDs{7},'_axon'])],col(1,:),allPreSynapse{7}, allPostSynapse{7});
DisplayTree(allTrees{21},[1],false, [eval([cellIDs{21},'_axon'])],col(2,:),allPreSynapse{21}, allPostSynapse{21});
set(gca,'color',[ctrans, 0.2]);
set(gcf,'color','none', 'units','normalized','outerposition',[0 0 1 1]);
h1 = gcf;
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

col = distinguishable_colors(numel(cellIDsDbx),cdbx);
DisplayTree(allTrees{2},[1],false, [eval([cellIDs{2},'_axon'])],col(1,:),allPreSynapse{2}, allPostSynapse{2});
DisplayTree(allTrees{9},[1],false, [eval([cellIDs{9},'_axon'])],col(4,:),allPreSynapse{9}, allPostSynapse{9});
set(gca,'color',[cdbx, 0.2]);
set(gcf,'color','none', 'units','normalized','outerposition',[0 0 1 1]);
h1 = gcf;
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

col = distinguishable_colors(numel(cellIDsL),cbarhl);
DisplayTree(allTrees{1},[1],false, [eval([cellIDs{1},'_axon'])],col(1,:),allPreSynapse{1}, allPostSynapse{1});
DisplayTree(allTrees{20},[1],false, [eval([cellIDs{20},'_axon'])],col(6,:),allPreSynapse{20}, allPostSynapse{20});
set(gca,'color',[cbarhl, 0.2]);
set(gcf,'color','none', 'units','normalized','outerposition',[0 0 1 1]);
h1 = gcf;
axis vis3d;
[h2,h3] = PlotViews(h1);
close(figure(4));

export_fig(h1,sprintf('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/ZfishpaperFigures/NewFigure4/AllBarhlXY.png'),'-r300','-transparent');
export_fig(h2,sprintf('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/ZfishpaperFigures/NewFigure4/AllBarhlYZ.png'),'-r300','-transparent');
export_fig(h3,sprintf('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/ZfishpaperFigures/NewFigure4/AllBarhlXZ.png'),'-r300','-transparent');

pause(100);
close all;
clear col;

