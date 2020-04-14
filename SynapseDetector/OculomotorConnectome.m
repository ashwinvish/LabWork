%%
load connMatCKT.mat; % Ordered connectivity matrix , rows are dendrites and columns are axons
load connMatOrderedTotalPostSynapses.mat; % vector with total number of input synapses along the rows
load groupEnds.mat % boundaries for the different groups
load groupIDs.mat % group IDs for neurons

% elements are total number of synapses




% assemble the ordered Matrix
[m,n,v] = intersect(mlOrdered,AllCells,'stable');

for i =1:size(m,1)
    mlConn(i,:) = ConnMatrixPre(v(i),v);
end   



subplot(2,2,1);
%heatmap(connMat,'GridVisible','off','Colormap',colorcet('L4'));
cspy(connMatCKT,'Colormap',colorcet('R3','N',15),'Levels',15,'MarkerSize',7);
axis square;
% colormap(colorcet('L17','N',15,'reverse',1));
lighting phong;
material shiny;
set(gca, 'Xtick',[],'XTickLabel',[],'YTick',[],'YTickLabel',[]);
set(gca,'color',[0.90,0.9,0.90]);
colorbar;
box on;


hold on
for i = 1:size(groupEnds,2)-1
line([0,size(connMatCKT,1)],[groupEnds(i),groupEnds(i)],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([groupEnds(i),groupEnds(i)],[0,size(connMatCKT,1)],'color',[0.5,0.5,0.5],'lineWidth',0.5);
end



subplot(2,2,2);
connMatNorm = bsxfun(@rdivide,connMatCKT,connMatOrderedTotalPostSynapses);
connMatNorm(connMatNorm>0.05) = 0.05;
cspy(connMatNorm ,'Colormap',colorcet('R3','N',15),'Levels',15,'MarkerSize',7);
caxis([0,0.05]);
axis square;
%colormap(colorcet('L17','N',15,'reverse',1));
lighting phong;
material shiny;
set(gca, 'Xtick',[],'XTickLabel',[],'YTick',[],'YTickLabel',[]);
set(gca,'color',[0.9,0.9,0.9]);
colorbar;
box on;

hold on
for i = 1:size(groupEnds,2)-1
line([0,size(connMatCKT,1)],[groupEnds(i),groupEnds(i)],'color',[0.5,0.5,0.5],'lineWidth',0.5);
line([groupEnds(i),groupEnds(i)],[0,size(connMatCKT,1)],'color',[0.5,0.5,0.5],'lineWidth',0.5);
end


